#' visual_magnitude
#' @description
#' Visual Magnitude quantifies the extent of a visible region
#' as perceived by an observer. It is derived from the surface's slope and
#' angle features, alongside the observer's relative distance from the area
#' (Chamberlain & Meitner).
#'
#' @param viewshed Viewshed object.
#' @param dsm SpatRaster. A digital surface / elevation model
#'
#' @return SpatRaster
#'
#' @examples
#' \donttest{
#' # Load a viewpoint
#' test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
#' # load dsm raster
#' dsm <- terra::rast(system.file("test_dsm.tif", package ="viewscape"))
#' # Compute viewshed
#' viewshed <- viewscape::compute_viewshed(dsm = dsm,
#'                                         viewpoints = test_viewpoint,
#'                                         offset_viewpoint = 6)
#' # Compute visual magnitude
#' vm <- viewscape::visual_magnitude(viewshed, dsm)
#' }
#'
#' @importFrom terra terrain
#'
#' @seealso [compute_viewshed()]
#'
#' @references Chamberlain, B. C., & Meitner, M. J. (2013).
#' A route-based visibility analysis for landscape management.
#' Landscape and Urban Planning, 111, 13-24.
#'
#' @export

visual_magnitude <- function(viewshed, dsm) {
  if (missing(dsm)) {
    stop("dsm is missing")
  }
  if (missing(viewshed)) {
    stop("Viewshed object is missing")
  }
  v <- viewshed
  # crop raster
  subdsm <- terra::crop(dsm, terra::ext(viewshed@extent, xy=TRUE))
  # get slope
  slope <- terra::terrain(subdsm, v="slope", neighbors=8, unit="radians")
  aspect <- terra::terrain(subdsm, v="aspect", neighbors=8, unit="radians")
  # direction <- terra::terrain(subdsm, v="flowdir")
  # convert raster to matrix
  dsm_matrix <- terra::as.matrix(subdsm, wide=TRUE)
  slope_matrix <- terra::as.matrix(slope, wide=TRUE)
  slope_matrix[is.nan(slope_matrix)] = 0
  aspect_matrix <- terra::as.matrix(aspect, wide=TRUE)
  aspect_matrix[is.nan(aspect_matrix)] = 0
  # compute visual magnitude
  vm_matrix <- VM(viewshed@visible, dsm_matrix,
                  slope_matrix, aspect_matrix,
                  viewshed@viewpos, viewshed@viewpoint[3],
                  viewshed@resolution[1])
  # v@visible <- vm_matrix
  # vm raster
  vm <- terra::rast(vm_matrix)
  # vmpoints <- filter_invisible(v, FALSE)
  # m <- terra::vect(sp::SpatialPoints(vmpoints))
  terra::crs(vm) <- viewshed@crs
  terra::ext(vm) <- terra::ext(viewshed@extent, xy = TRUE)
  vm <- terra::classify(vm, cbind(-9, NA))
  # vm <- terra::mask(filter_invisible(v, TRUE), m)
  return(vm)
}
