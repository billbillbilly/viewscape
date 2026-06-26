#' pano_view
#'
#' @description
#' Generate an equirectangular or cylindrical panoramic view from a DSM.
#'
#' For `method = "equirectangular"`:
#' \itemize{
#'   \item \strong{No semantic}: automatically downloads ESRI WorldImagery tiles
#'     via \pkg{greenSD} and projects the satellite RGB colours onto the
#'     panoramic rays, producing a 3-band colour panorama (R, G, B, 0–255).
#'   \item \strong{With semantic}: projects the supplied land-use / land-cover
#'     raster onto the panoramic rays, returning a 2-layer SpatRaster (depth
#'     + land-cover codes).
#' }
#' For `method = "cylindrical"` the function always returns a single-layer
#' distance panorama (with an optional second semantic layer).
#'
#' @param dsm A SpatRaster object (single layer) of surface elevation.
#' @param vpt A POINT sf object or numeric coordinate pair representing the
#'   viewpoint (in projected coordinates).
#' @param h Numeric. Height of the viewpoint above the ground surface (in
#'   metres). Default: 6.
#' @param semantic A SpatRaster object (single layer) of land-use / land-cover
#'   class values. When supplied with `method = "equirectangular"`, the
#'   panorama encodes the land-cover code at each ray-hit cell (depth layer
#'   plus semantic layer). Binary 0/1 rasters are treated as occupied vertical
#'   structures. When supplied with `method = "cylindrical"`, the return value
#'   has two layers: depth and semantic class.
#' @param method Character. `"equirectangular"` (default) or `"cylindrical"`.
#' @param sky_threshold Numeric. Elevation buffer (metres) below the ray
#'   position that is still treated as terrain (default: 3.0). Must be less
#'   than \code{h}; if not, it is automatically clamped to \code{h - 0.01}.
#' @param step_size Numeric. Ray marching step size in DSM coordinate units
#'   (default: 0.5).
#' @param max_dist Numeric. Maximum ray distance (default: 500).
#' @param pano_dim Integer vector of length 2: panorama height and width in
#'   pixels (default: `c(256, 512)`).
#' @param ares Numeric. Cylindrical vertical angular resolution in radians per
#'   pixel. Ignored when `method = "equirectangular"` (default: 0.0174533).
#' @param heading Numeric. Camera heading in degrees (default: 0 = north).
#'   See Details.
#' @param color_source Optional SpatRaster with exactly 3 bands (R, G, B,
#'   values 0–255). Overrides the automatic greenSD satellite download when
#'   `method = "equirectangular"` and `semantic = NULL`. Useful when you
#'   already have a local satellite image.
#' @param sky_color Numeric vector of length 3 (R, G, B, 0–255). Colour
#'   painted on sky pixels in the colour panorama. Default: sky blue
#'   `c(135, 206, 235)`.
#' @param satellite_zoom Integer. Tile zoom level passed to
#'   `greenSD::get_tile_green()` when auto-fetching satellite imagery
#'   (default: 17, approximately 1.2 m/px at 42 degrees N).
#' @param plot Logical. Whether to plot the result (default: `FALSE`).
#' @param legend Logical. Whether to display a legend when plotting
#'   (default: `TRUE`).
#' @param axes Logical. Whether to display axes when plotting (default: `TRUE`).
#'
#' @return
#' \itemize{
#'   \item `method = "equirectangular"`, no `semantic`: 3-band SpatRaster
#'     (R, G, B, 0–255). Plot with `terra::plotRGB()`.
#'   \item `method = "equirectangular"`, with `semantic`: 2-layer SpatRaster
#'     (depth in metres, land-cover code).
#'   \item `method = "cylindrical"`, no `semantic`: single-layer distance
#'     SpatRaster (`NA` = sky).
#'   \item `method = "cylindrical"`, with `semantic`: 2-layer SpatRaster
#'     (depth, semantic class).
#' }
#'
#' @details
#' For `heading`:
#' \itemize{
#'   \item `heading = 0`   → facing north (default)
#'   \item `heading = 90`  → facing east
#'   \item `heading = 180` → facing south
#'   \item `heading = 270` → facing west
#' }
#'
#' The equirectangular colour panorama requires the \pkg{greenSD} package
#' (install with `devtools::install_github("billbillbilly/greenSD")`).
#' It calls `greenSD::get_tile_green(bbox, zoom, provider = "esri")` and
#' uses the raw RGB tile (`[[2]]`) as the colour source.
#'
#' @examples
#' dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
#' vpt <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
#' # cylindrical distance panorama (no internet needed)
#' result <- viewscape::pano_view(dsm, vpt, h = 6, method = "cylindrical")
#'
#' @importFrom terra colFromX rowFromY as.matrix rast extract plot ext resample
#'   global nlyr plotRGB xmin xmax ymin ymax crs project
#' @importFrom sf st_coordinates
#' @importFrom grDevices gray.colors
#'
#' @references
#' Lu, X., Li, Z., Cui, Z., Oswald, M. R., Pollefeys, M., & Qin, R. (2020).
#' Geometry-aware satellite-to-ground image synthesis for urban areas.
#' In Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern
#' Recognition (pp. 859-867).
#' @export
pano_view <- function(dsm = NULL,
                      vpt = NULL,
                      h = 6,
                      semantic = NULL,
                      method = 'equirectangular',
                      sky_threshold = 3.0,
                      step_size = 0.5,
                      max_dist = 500,
                      pano_dim = c(256, 512),
                      ares = 0.0174533,
                      heading = 0,
                      color_source = NULL,
                      sky_color = c(135, 206, 235),
                      satellite_zoom = 17L,
                      plot = FALSE,
                      legend = TRUE,
                      axes = TRUE) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (is.null(dsm)) stop("dsm is missing")
  if (is.null(vpt)) stop("vpt is missing")
  method <- match.arg(method, c("cylindrical", "equirectangular"))
  if (!is.numeric(h) || length(h) != 1 || is.na(h))
    stop("h must be a single numeric value")
  if (!is.numeric(pano_dim) || length(pano_dim) != 2 ||
      any(is.na(pano_dim)) || any(pano_dim <= 0))
    stop("pano_dim must be a positive numeric vector of length 2")
  if (!is.numeric(step_size) || length(step_size) != 1 ||
      is.na(step_size) || step_size <= 0)
    stop("step_size must be a single positive numeric value")
  if (!is.numeric(max_dist) || length(max_dist) != 1 ||
      is.na(max_dist) || max_dist <= 0)
    stop("max_dist must be a single positive numeric value")
  if (!is.numeric(ares) || length(ares) != 1 || is.na(ares) || ares <= 0)
    stop("ares must be a single positive numeric value")
  if (!is.numeric(sky_color) || length(sky_color) != 3)
    stop("sky_color must be a numeric vector of length 3 (R, G, B)")
  sky_color <- pmin(pmax(sky_color, 0), 255)

  # ── Validate semantic ───────────────────────────────────────────────────────
  if (!is.null(semantic)) {
    if (!inherits(semantic, "SpatRaster"))
      stop("semantic must be a SpatRaster")
    if (dim(semantic)[3] != 1)
      stop("semantic must have exactly one layer")
    if (!terra::hasValues(semantic))
      stop("semantic has no values")
    if (!terra::compareGeom(dsm, semantic, stopOnError = FALSE)) {
      warning("resampling semantic to the dsm grid",
              call. = FALSE)
      semantic <- terra::resample(semantic, dsm, method = "near")
    }
    if (!identical(dim(dsm)[1:2], dim(semantic)[1:2]))
      stop("semantic must resolve to the same dimensions as dsm")
    if (all(is.na(terra::values(semantic))))
      stop("semantic has no non-NA values after alignment to dsm")
  }

  # ── Validate explicit color_source ─────────────────────────────────────────
  if (!is.null(color_source)) {
    if (!inherits(color_source, "SpatRaster"))
      stop("color_source must be a SpatRaster")
    if (terra::nlyr(color_source) != 3)
      stop("color_source must have exactly 3 bands (R, G, B)")
    if (!terra::compareGeom(dsm, color_source, stopOnError = FALSE))
      color_source <- terra::resample(color_source, dsm, method = "bilinear")
  }

  # ── Viewpoint coordinates and elevation ────────────────────────────────────
  if (is.numeric(vpt)) {
    if (length(vpt) < 2) stop("vpt must contain x and y coordinates")
    vpt_coords <- vpt[1:2]
  } else {
    vpt_coords <- sf::st_coordinates(vpt)[1, ]
  }
  vpt_row <- terra::rowFromY(dsm, vpt_coords[2])
  vpt_col <- terra::colFromX(dsm, vpt_coords[1])
  .ex <- terra::extract(dsm, matrix(vpt_coords, ncol = 2))
  terrain_val <- .ex[1, ncol(.ex)]
  vpt_z <- terrain_val + h
  if (is.na(vpt_row) || is.na(vpt_col) || is.na(vpt_z))
    stop("vpt must be located on a non-NA cell of dsm")

  # Clamp sky_threshold so it is strictly below h.
  # If sky_threshold >= h, every ray hits the ground at the very first step
  # (self-occlusion) because the observer is already inside the threshold band.
  if (sky_threshold >= h) {
    sky_threshold <- max(0.01, h - 0.01)
    message("sky_threshold clamped to ", round(sky_threshold, 3),
            " (must be < h = ", h, ")")
  }

  dsm_matrix      <- terra::as.matrix(dsm, wide = TRUE)
  semantic_matrix <- if (is.null(semantic)) NULL else
                       terra::as.matrix(semantic, wide = TRUE)
  semantic_is_binary_mask <- FALSE
  if (!is.null(semantic)) {
    sr <- terra::global(semantic, c("min", "max"), na.rm = TRUE)
    semantic_is_binary_mask <- isTRUE(sr[1, "min"] >= 0 && sr[1, "max"] <= 1)
  }
  dsm_resolution  <- terra::res(dsm)
  orientation_rad <- heading * pi / 180 + pi / 2

  # ── Helper: run color raycaster ────────────────────────────────────────────
  .color_pano <- function(sat) {
    sat_r <- terra::as.matrix(sat[[1]], wide = TRUE)
    sat_g <- terra::as.matrix(sat[[2]], wide = TRUE)
    sat_b <- terra::as.matrix(sat[[3]], wide = TRUE)
    pano  <- dsm_to_color_equirectangular(
      dsm_matrix, sat_r, sat_g, sat_b,
      vpt_x         = vpt_col - 1,
      vpt_y         = vpt_row - 1,
      vpt_z         = vpt_z,
      pano_height   = as.integer(pano_dim[1]),
      pano_width    = as.integer(pano_dim[2]),
      xres          = dsm_resolution[1],
      yres          = dsm_resolution[2],
      step_size     = step_size,
      max_dist      = max_dist,
      sky_threshold = sky_threshold,
      orientation   = orientation_rad,
      sky_r         = sky_color[1],
      sky_g         = sky_color[2],
      sky_b         = sky_color[3]
    )
    r <- c(terra::rast(pano$r), terra::rast(pano$g), terra::rast(pano$b))
    names(r) <- c("R", "G", "B")
    terra::ext(r) <- terra::ext(0, 360, -90, 90)
    if (plot) terra::plotRGB(r, r = 1, g = 2, b = 3, axes = axes)
    r
  }

  # ── Dispatch ────────────────────────────────────────────────────────────────
  if (method == "cylindrical") {

    if (is.null(semantic_matrix)) {
      pano <- dsm_to_pano(dsm_matrix,
                          vpt_x         = vpt_col - 1,
                          vpt_y         = vpt_row - 1,
                          vpt_z         = vpt_z,
                          pano_height   = as.integer(pano_dim[1]),
                          pano_width    = as.integer(pano_dim[2]),
                          xres          = dsm_resolution[1],
                          yres          = dsm_resolution[2],
                          ares          = ares,
                          step_size     = step_size,
                          max_dist      = max_dist,
                          sky_value     = NA_real_,
                          sky_threshold = sky_threshold,
                          orientation   = orientation_rad)
    } else {
      pano <- dsm_semantic_to_pano(dsm_matrix,
                                   semantic_matrix,
                                   vpt_x         = vpt_col - 1,
                                   vpt_y         = vpt_row - 1,
                                   vpt_z         = vpt_z,
                                   pano_height   = as.integer(pano_dim[1]),
                                   pano_width    = as.integer(pano_dim[2]),
                                   xres          = dsm_resolution[1],
                                   yres          = dsm_resolution[2],
                                   ares          = ares,
                                   step_size     = step_size,
                                   max_dist      = max_dist,
                                   sky_value     = NA_real_,
                                   sky_threshold = sky_threshold,
                                   orientation   = orientation_rad)
    }

  } else {  # equirectangular ─────────────────────────────────────────────────

    if (is.null(semantic_matrix)) {
      # ── Color panorama (satellite RGB) ──────────────────────────────────────
      if (!is.null(color_source)) {
        # Explicit override supplied by user
        return(.color_pano(color_source))
      }

      # Auto-fetch ESRI WorldImagery tiles via greenSD
      if (!requireNamespace("greenSD", quietly = TRUE)) {
        stop("Package 'greenSD' is required for equirectangular colour panoramas.\n",
             "Install with: devtools::install_github(\"billbillbilly/greenSD\")\n",
             "Or supply a pre-downloaded raster via the `color_source` argument.")
      }
      message("Fetching ESRI WorldImagery tiles via greenSD (zoom = ",
              satellite_zoom, ")...")
      dsm_geo <- terra::project(dsm, "EPSG:4326")
      bbox    <- c(terra::xmin(dsm_geo), terra::ymin(dsm_geo),
                   terra::xmax(dsm_geo), terra::ymax(dsm_geo))
      tiles   <- greenSD::get_tile_green(bbox     = bbox,
                                         zoom     = satellite_zoom,
                                         provider = "esri")
      sat_rgb <- tiles[[2]]                                    # raw RGB raster
      sat_rgb <- terra::project(sat_rgb, terra::crs(dsm))
      sat_rgb <- terra::resample(sat_rgb, dsm, method = "bilinear")
      return(.color_pano(sat_rgb))

    } else {
      # ── Semantic / land-cover panorama ──────────────────────────────────────
      pano <- dsm_semantic_to_equirectangular(
        dsm_matrix,
        semantic_matrix,
        vpt_x            = vpt_col - 1,
        vpt_y            = vpt_row - 1,
        vpt_z            = vpt_z,
        ground_z         = terrain_val,
        pano_height      = as.integer(pano_dim[1]),
        pano_width       = as.integer(pano_dim[2]),
        xres             = dsm_resolution[1],
        yres             = dsm_resolution[2],
        step_size        = step_size,
        max_dist         = max_dist,
        sky_value        = NA_real_,
        sky_threshold    = sky_threshold,
        orientation      = orientation_rad,
        extrude_semantic = semantic_is_binary_mask
      )
    }

  }

  # ── Assemble non-colour output ───────────────────────────────────────────────
  if (is.null(semantic_matrix)) {
    r <- terra::rast(pano)
    names(r) <- paste0(method, "_depth")
  } else {
    r <- c(terra::rast(pano$depth), terra::rast(pano$semantic))
    names(r) <- c(paste0(method, "_depth"), paste0(method, "_semantic"))
  }
  terra::ext(r) <- terra::ext(0, 360, -90, 90)

  if (plot) {
    if (is.null(semantic_matrix)) {
      terra::plot(r,
                  col    = gray.colors(100, start = 0, end = 1),
                  legend = legend,
                  axes   = axes)
    } else {
      old_par <- graphics::par(mfrow = c(1, 2))
      on.exit(graphics::par(old_par), add = TRUE)
      terra::plot(r[[1]],
                  col    = gray.colors(100, start = 0, end = 1),
                  legend = legend, axes = axes,
                  main   = names(r)[1])
      if (semantic_is_binary_mask) {
        terra::plot(r[[2]],
                    col    = c("gray50", "black"),
                    breaks = c(-Inf, 0.5, Inf),
                    colNA  = "white",
                    legend = legend, axes = axes,
                    main   = names(r)[2])
      } else {
        terra::plot(r[[2]], legend = legend, axes = axes, main = names(r)[2])
      }
    }
  }
  return(r)
}
