#' calculate_feature
#' @description The calculate_feature function is designed to extract specific
#' feature-related information within a viewshed. It allows you to compute
#' the proportion of the feature that is present in the viewshed.
#'
#' @param viewshed Viewshed object.
#' @param feature Raster. Land cover or land use
#' @param type Numeric. The input type of land cover raster.
#' type=1: percentage raster (that represents the percentage of
#' area in each cell).
#' type=2: binary raster (that only uses two values to represent whether
#' the feature exists in each cell).
#' @param exclude_value Numeric. the value of those cells need to be excluded
#' in the analysis. If type = 2, exclude_value is reqired.
#'
#' @return Numeric. The canopy area in the viewshed.
#' @export
#'
#' @examples
#' library(viewscape)
#' # Load a viewpoint
#' test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
#' # load dsm raster
#' dsm <- terra::rast(system.file("test_dsm.tif", package ="viewscape"))
#' #Compute viewshed
#' viewshed <- viewscape::compute_viewshed(dsm = dsm,
#'                                         viewpoints = test_viewpoint,
#'                                         offset_viewpoint = 6)
#' # load canopy raster
#' test_canopy <- terra::rast(system.file("test_canopy.tif",
#'                                        package ="viewscape"))
#' # calculate the percentage of canopy coverage
#' test_canopy_proportion <- viewscape::calculate_feature(viewshed = viewshed,
#'                                                        feature = test_canopy,
#'                                                        type = 2,
#'                                                        exclude_value = 0)

calculate_feature <- function(viewshed,
                              feature,
                              type,
                              exclude_value){
  if (missing(viewshed)) stop("viewshed is missing")
  if (missing(feature))  stop("feature is missing")
  if (missing(type))     stop("type is missing")
  if (type == 2 && missing(exclude_value)) stop("please specify exclude_value")

  if (isFALSE(terra::crs(feature, proj = TRUE) == viewshed@crs)) {
    feature <- terra::project(feature, y = viewshed@crs)
  }

  pt <- filter_invisible(viewshed, FALSE)
  feature_df <- terra::extract(feature, pt)
  feature_df <- feature_df[, ncol(feature_df)]

  if (type == 1) {
    output <- sum(feature_df, na.rm = TRUE) / length(pt[, 1])
  } else if (type == 2) {
    feature_ <- feature_df[feature_df != exclude_value]
    output <- length(feature_) / length(pt[, 1])
  }
  return(output)
}
