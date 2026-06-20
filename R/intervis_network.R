#' intervis_network
#' @description Computes an intervisibility network among a set of viewpoints,
#' determining which pairs of viewpoints have a direct line of sight to each other.
#' Results can be returned as a binary adjacency matrix or as an sf collection of
#' lines suitable for mapping the network.
#'
#' @param viewpoints sf point(s), or a matrix/data.frame of x,y coordinates.
#'   At least 2 viewpoints are required.
#' @param dsm SpatRaster. A digital surface / elevation model.
#' @param offset_viewpoint Numeric. Observer height above the surface **in metres**.
#'   Either a single value applied to all viewpoints, or a vector of
#'   per-viewpoint heights. Default is 1.7 m. Automatically converted to the
#'   CRS unit of the DSM.
#' @param r Numeric (optional). Analysis radius **in metres**. Pairs separated
#'   by more than \code{r} are treated as non-visible. Default is 1000 m.
#'   Automatically converted to the CRS unit of the DSM.
#' @param refraction_factor Numeric. Atmospheric refraction coefficient.
#'   Default is 0.13.
#' @param method Character. Viewshed algorithm: \code{"plane"} (default) or
#'   \code{"los"}. See \code{\link{compute_viewshed}} for details.
#' @param output Character. Return type: \code{"matrix"} (default) returns an
#'   N x N integer adjacency matrix where entry [i,j] is 1 if viewpoint i can
#'   see viewpoint j, and 0 otherwise (NA on the diagonal); \code{"lines"}
#'   returns an sf data frame of LINESTRINGs connecting visible pairs, with
#'   columns \code{from}, \code{to}, and \code{mutual} (TRUE when both
#'   directions are visible).
#'
#' @return An N x N integer matrix (\code{output = "matrix"}) or an sf data
#'   frame of lines (\code{output = "lines"}).
#'
#' @examples
#' \donttest{
#' library(viewscape)
#' # Load viewpoints
#' test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
#' # Load DSM
#' dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
#' # Compute intervisibility matrix
#' mat <- intervis_network(test_viewpoint, dsm, offset_viewpoint = 1.7)
#' # Compute intervisibility as lines
#' net <- intervis_network(test_viewpoint, dsm, offset_viewpoint = 1.7,
#'                         output = "lines")
#' }
#'
#' @seealso [compute_viewshed()] [fov_mask()]
#' @export

intervis_network <- function(viewpoints,
                              dsm,
                              offset_viewpoint = 1.7,
                              r = NULL,
                              refraction_factor = 0.13,
                              method = "plane",
                              output = "matrix") {
  # --- input validation ---
  if (missing(viewpoints)) stop("viewpoints is missing")
  if (missing(dsm))        stop("dsm is missing")
  if (!output %in% c("matrix", "lines"))
    stop("output must be 'matrix' or 'lines'")
  if (!method %in% c("plane", "los"))
    stop("method must be 'plane' or 'los'")

  if (inherits(viewpoints, "sf")) {
    vp_coords <- sf::st_coordinates(viewpoints)[, 1:2, drop = FALSE]
  } else if (is.matrix(viewpoints) || is.data.frame(viewpoints)) {
    vp_coords <- as.matrix(viewpoints)[, 1:2, drop = FALSE]
  } else {
    stop("viewpoints must be an sf object or a matrix/data.frame of x,y coordinates")
  }

  n <- nrow(vp_coords)
  if (n < 2) stop("At least 2 viewpoints are required")

  if (length(offset_viewpoint) == 1) {
    offset_viewpoint <- rep(offset_viewpoint, n)
  } else if (length(offset_viewpoint) != n) {
    stop("offset_viewpoint must be length 1 or the same length as the number of viewpoints")
  }

  # Convert metre-based inputs to the DSM's native CRS unit (e.g. feet)
  if (is.null(r)) r <- 1000
  r                <- m_to_crs_units(r,    dsm)
  offset_viewpoint <- m_to_crs_units(offset_viewpoint, dsm)

  # --- N x N adjacency matrix (NA on diagonal) ---
  result <- matrix(NA_integer_, nrow = n, ncol = n)

  for (i in seq_len(n)) {
    vp_i <- c(vp_coords[i, 1], vp_coords[i, 2])
    vs_i <- radius_viewshed(dsm, r, refraction_factor, vp_i,
                            offset_viewpoint[i], 0, method)
    vs_rast <- filter_invisible(vs_i, TRUE)

    for (j in seq_len(n)) {
      if (i == j) next
      target_pt <- matrix(c(vp_coords[j, 1], vp_coords[j, 2]), nrow = 1)
      ex  <- terra::extract(vs_rast, target_pt)
      val <- ex[1, ncol(ex)]
      result[i, j] <- if (!is.na(val) && val > 0) 1L else 0L
    }
  }

  if (output == "matrix") {
    return(result)
  }

  # --- build sf LINESTRING collection for visible pairs ---
  lines_list <- vector("list", 0)
  from_ids   <- integer(0)
  to_ids     <- integer(0)
  mutual     <- logical(0)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (j <= i) next  # upper triangle only — avoids duplicate lines
      vis_ij <- isTRUE(result[i, j] == 1L)
      vis_ji <- isTRUE(result[j, i] == 1L)
      if (vis_ij || vis_ji) {
        lines_list <- c(lines_list,
                        list(sf::st_linestring(rbind(vp_coords[i, ],
                                                     vp_coords[j, ]))))
        from_ids <- c(from_ids, i)
        to_ids   <- c(to_ids,   j)
        mutual   <- c(mutual, vis_ij && vis_ji)
      }
    }
  }

  if (length(lines_list) == 0) {
    warning("No intervisible pairs found")
    return(sf::st_sf(from    = integer(0),
                     to      = integer(0),
                     mutual  = logical(0),
                     geometry = sf::st_sfc(crs = sf::st_crs(dsm))))
  }

  sf::st_sf(
    from     = from_ids,
    to       = to_ids,
    mutual   = mutual,
    geometry = sf::st_sfc(lines_list, crs = sf::st_crs(dsm))
  )
}
