#' calculate_diversity
#' @description The calculate_diversity function is designed to calculate landscape
#' diversity metrics within a viewshed. It takes as input a land cover raster,
#' a viewshed object representing the observer's line of sight, and an optional
#' parameter to compute class proportions.
#'
#' @param viewshed Viewshed object.
#' @param land Raster. The raster of land use/land cover representing different
#' land use/cover classes.
#' @param proportion logical (Optional), indicating whether to return class
#' proportions along with the Shannon Diversity Index (SDI). (default is FALSE).
#' @return List. a list containing the Shannon Diversity Index (SDI) and,
#' if the proportion parameter is set to TRUE, a table of class proportions
#' within the viewshed.
#' @import terra
#' @importFrom dplyr count
#' @importFrom dplyr select
#' @importFrom dplyr %>%
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' library(viewscape)
#' # Load a viewpoint
#' test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
#' # load dsm raster
#' dsm <- terra::rast(system.file("test_dsm.tif", package ="viewscape"))
#' #Compute viewshed
#' output <- viewscape::compute_viewshed(dsm = dsm,
#'                                       viewpoints = test_viewpoint,
#'                                       offset_viewpoint = 6, r = 1600)
#' # load landuse raster
#' test_landuse <- terra::rast(system.file("test_landuse.tif",
#'                                         package ="viewscape"))
#' diversity <- viewscape::calculate_diversity(output,
#'                                             test_landuse)
#'

calculate_diversity <- function(viewshed,
                                land,
                                proportion = FALSE){
  if (isFALSE(terra::crs(land, proj = TRUE) == viewshed@crs)) {
    land <- terra::project(land, y = viewshed@crs)
  }
  pt <- filter_invisible(viewshed, FALSE)
  land <- terra::crop(land, terra::ext(viewshed@extent, xy = TRUE))
  # calculate the proportion of each class
  land_class <- terra::extract(land, pt)
  land_class <- land_class[, ncol(land_class)]
  land_class <- as.data.frame(land_class)
  colnames(land_class)[1] <- "type"
  land_class <- land_class[!is.na(land_class$type), , drop = FALSE]
  land_class <- dplyr::count(land_class, .data$type)
  # land_class <- dplyr::count(land_class, "type")
  total <- sum(land_class$n)
  land_class$proportion <- land_class$n/total
  # calculate Shannon diversity index
  p <- land_class$proportion
  sdi <- sd_index(p)
  if (proportion) {
    sub_land_class <- land_class %>%
      dplyr::select(.data$type, proportion)
    return(list(SDI=sdi, Proportion=sub_land_class))
  } else {
    return(sdi)
  }
}
