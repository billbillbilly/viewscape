#' pano_view
#'
#' @description
#' Generate tangential panoramic view from a DSM based on the Geo-transformation
#' method proposed by Lu et al. (2020).
#'
#' @param dsm A SpatRaster object (single layer) of surface elevation.
#' @param semantic A SpatRaster object (single layer) of semantic values such as
#' land cover and land use.
#' @param vpt A POINT sf object representing the viewpoint
#' (in projected coordinates).
#' @param h Numeric. Height of the viewpoint above the ground surface (in meters).
#' @param method character. `cylindrical` or `equirectangular`
#' @param sky_threshold Numeric. Elevation buffer (in meters) to consider a ray
#' reaching sky (default: 3.0).
#' @param step_size Numeric. Distance (in pixels) between sampling steps
#' (default: 0.5).
#' @param max_dist Numeric. Maximum ray distance (default: 500).
#' @param pano_dim Integer vector of size 2: height and width of the output
#' panorama (default: c(256, 512)).
#' @param heading numeric. Horizontal direction in degree. (default: 0). See Details.
#' @param plot logical. Whether to plot result. (default: FLASE)
#' @param legend logical. Whether to display lengend when plotting the result.
#' (default: TRUE)
#' @param axes logical. Whether to display axis when plotting the result.
#' (default: TRUE)
#'
#' @return A SpatRaster object representing the distance panorama, with `NA`
#' where sky is detected.
#'
#' @details
#' For `heading`:
#' \itemize{
#'  \item heading = 0 → facing north (default)
#'  \item heading = 90 → facing east
#'  \item heading = 180 → facing south
#'  \item heading = 270 → facing due west
#' }
#'
#' @examples
#' dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
#' vpt <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
#' result <- streetscape::pano_view(dsm, vpt, 6, method = 'cylindrical')
#'
#' @importFrom terra colFromX rowFromY as.matrix rast extract vect plot values
#' @importFrom sf st_coordinates
#' @importFrom torch torch_tensor torch_float as_array
#'
#' @references
#' Lu, X., Li, Z., Cui, Z., Oswald, M. R., Pollefeys, M., & Qin, R. (2020).
#' Geometry-aware satellite-to-ground image synthesis for urban areas.
#' In Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern
#' Recognition (pp. 859-867).
#' @export
pano_view <- function(dsm = NULL,
                      semantic = NULL,
                      vpt = NULL,
                      h = 6,
                      method = 'equirectangular',
                      sky_threshold = 3.0,
                      step_size = 0.5,
                      max_dist = 500,
                      pano_dim = c(256, 512),
                      heading = 0,
                      plot = FALSE,
                      legend = TRUE,
                      axes = TRUE) {

  # Extract coordinates and elevation at the viewpoint
  vpt_coords <- sf::st_coordinates(vpt)[1, ]
  vpt_row <- terra::rowFromY(dsm, vpt_coords[2])
  vpt_col <- terra::colFromX(dsm, vpt_coords[1])
  terrain_val <- terra::extract(dsm, matrix(vpt_coords, ncol = 2))[[1]]
  vpt_z <- terrain_val + h
  if (method == 'cylindrical') {
    dsm_matrix <- terra::as.matrix(dsm, wide = TRUE)
    # compute orientation
    orientation_rad <- heading * pi / 180 + pi/2

    pano <- dsm_to_pano(dsm_matrix,
                        vpt_x = vpt_col,
                        vpt_y = vpt_row,
                        vpt_z = vpt_z,
                        pano_height = pano_dim[1],
                        pano_width = pano_dim[2],
                        step_size = step_size,
                        max_dist = max_dist,
                        sky_value = NA_real_,
                        sky_threshold = sky_threshold,
                        orientation = orientation_rad)

    # Convert to SpatRaster with dummy extents
    r <- terra::rast(pano)
    r <- r * step_size
  } else if (method == 'equirectangular') {
    dsm_tensor <- torch::torch_tensor(array(terra::values(dsm),
                                            dim = c(1, 1, nrow(dsm), ncol(dsm))
                                            )
                                      )$to(dtype = torch::torch_float())
    if (!is.null(semantic)) {
      sem_tensor <- torch::torch_tensor(array(terra::values(semantic),
                                              dim = c(1, 1,
                                                      nrow(semantic),
                                                      ncol(semantic)
                                                     )
                                              )
                                        )$to(dtype = torch::torch_float())
    }

    pano <- geo_projection(sate_depth = dsm_tensor,
                           sate_sem = if (is.null(semantic)) semantic else sem_tensor,
                           orientation = torch::torch_tensor(pi),
                           gsd = 1,
                           pano_size = pano_dim,
                           view_row = vpt_row,
                           view_col = vpt_col,
                           view_height = vpt_z)
    pano_depth <- torch::as_array(pano$pano_depth$squeeze())
    pano_sem <- torch::as_array(pano$pano_sem$squeeze())

    # Convert to terra raster
    r1 <- terra::rast(pano_depth)
    r2 <- terra::rast(pano_sem)
    r <- list(r1, r2)
  }

  # names(r) <- "tangential_distance"
  if (plot) {
    terra::plot(r, col = gray.colors(100, start = 0, end = 1),
                legend = legend,
                axes = axes)
  }
  return(r)
}

#### testing
# Load DSM and DTM from your 'viewscape' test files
dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
sem <- terra::rast(system.file("test_landuse.tif", package = "viewscape"))

# Load a test observer point
vpt <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
h <- 6
vpt_coords <- sf::st_coordinates(vpt)[1, ]
vpt_row <- terra::rowFromY(dsm, vpt_coords[2])
vpt_col <- terra::colFromX(dsm, vpt_coords[1])
vpt_z <- h

dsm_tensor <- torch::torch_tensor(array(terra::values(dsm - terra::minmax(dsm)[1]),
                                        dim = c(1, 1, nrow(dsm), ncol(dsm))
)
)$to(dtype = torch::torch_float())

sem_tensor <- torch::torch_tensor(array(terra::values(sem),
                                        dim = c(1, 1,
                                                nrow(sem),
                                                ncol(sem)
                                        )
)
)$to(dtype = torch::torch_float())

pano <- geo_projection(sate_depth = dsm_tensor,
                       sate_sem =  sem_tensor,
                       orientation = torch::torch_tensor(pi),
                       gsd = 0.1,
                       pano_size = c(256, 512),
                       view_row = vpt_row,
                       view_col = vpt_col,
                       view_height = vpt_z)
pano_depth <- torch::as_array(pano$pano_depth$squeeze())
pano_sem <- torch::as_array(pano$pano_sem$squeeze())

# Convert to terra raster
r1 <- terra::rast(pano_depth)
r2 <- terra::rast(pano_sem)
r <- list(r1, r2)
