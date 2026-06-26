#' compute_viewshed
#' @description The compute_viewshed function is designed for computing viewsheds,
#' which are areas visible from specific viewpoints, based on a Digital Surface
#' Model (DSM). It provides flexibility for single or multi-viewpoint analyses
#' and allows options for parallel processing, raster output, and plotting.
#'
#' @param dsm Raster, the digital surface model/digital elevation model
#' @param viewpoints sf point(s) or vector including x,y coordinates of a viewpoint
#' or a matrix including several viewpoints with x,y coordinates
#' @param offset_viewpoint numeric, observer height above the surface **in metres**
#' (default is 1.7 m). Automatically converted to the CRS unit of the DSM.
#' @param offset_height numeric, target height above the surface **in metres**
#' (default is 0 m). Automatically converted to the CRS unit of the DSM.
#' @param r Numeric (optional), analysis radius **in metres**
#' (default is 1000 m). Automatically converted to the CRS unit of the DSM.
#' @param refraction_factor Number, indicating the refraction factor.
#' The refraction factor adjusts the effect of atmospheric refraction
#' on the apparent curvature of the Earth. In most standard applications, a refraction factor
#' of 0.13 is used, and so does this function. However, the appropriate refraction factor may
#' vary depending on environmental conditions.
#' @param method Character, The algorithm for computing a viewshed:
#' "plane", "los", and "view_tree" (see details). "plane" is used as default.
#' @param parallel Logical, (default is FALSE) indicating if parallel computing
#' should be used to compute viewsheds of multiview points. When it is TRUE,
#' arguements 'raster' and 'plot' are ignored
#' @param workers Numeric, indicating the number of CPU cores that will be used
#' for parallel computing. It is required if 'parallel' is 'TRUE'.
#' @param raster Logical, (default is FALSE) if it is TRUE, the raster of
#' viewshed will be returned.
#' The default is FALSE
#' @param plot Logical, (default is FALSE) if it is TRUE, the raster of
#' viewshed will be displayed
#'
#' @return Raster or list. For single-viewpoint analysis, the function returns
#' either a raster (raster is TRUE) or a viewshed object. Value 1 means visible while
#' value 0 means invisible. For multi-viewpoint analysis, a list of viewsheds
#' is returned.
#'
#'
#’ @details For method, "plane" is the reference plane algorithm introduced by
#’ Wang et al. (2000) and "los" is the line of sight algorithm (Franklin & Ray, 1994).
#’ The reference plane algorithm can be more time-efficient than the line of sight
#’ algorithm, whereas the accuracy of the line of sight is better.
#’
#’ "view_tree" is the view-tree algorithm (Wang et al. 2022). It converts the DEM
#’ into a tree data structure (the view-tree) based on spatial occlusion relationships
#’ among cells, then traverses the tree depth-first, propagating a minimum visible
#’ slope (MVS) from parent to child nodes to determine visibility. The method achieves
#’ precision comparable to "los" while being approximately 40% faster than R2-class
#’ methods. Note that this is a pure-R implementation, so performance on very large
#’ DEMs will be slower than the compiled "plane" and "los" methods.
#’
#’ @seealso [fov_mask()] [visual_magnitude()]
#’
#’ @references Franklin, W. R., & Ray, C. (1994, May).
#’ Higher isn’t necessarily better: Visibility algorithms and experiments.
#’ In Advances in GIS research: sixth international symposium on spatial
#’ data handling (Vol. 2, pp. 751-770). Edinburgh: Taylor & Francis.
#’
#’ Wang, J., Robinson, G. J., & White, K. (2000).
#’ Generating viewsheds without using sightlines.
#’ Photogrammetric engineering and remote sensing, 66(1), 87-90.
#’
#’ Wang, Z., Xiong, L., Guo, Z., Zhang, W., & Tang, G. (2022).
#’ A view-tree method to compute viewsheds from digital elevation models.
#’ International Journal of Geographical Information Science.
#’ https://doi.org/10.1080/13658816.2022.2094385
#'
#' @useDynLib viewscape
#' @import pbmcapply
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @importFrom terra ext
#' @importFrom terra res
#' @importFrom terra plot
#' @importFrom parallel detectCores
#' @export
#'
#' @examples
#' # Load a viewpoint
#' test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
#' # load dsm raster
#' dsm <- terra::rast(system.file("test_dsm.tif", package ="viewscape"))
#' #Compute viewshed
#' output <- viewscape::compute_viewshed(dsm = dsm,
#'                                       viewpoints = test_viewpoint,
#'                                       offset_viewpoint = 6, r = 1600)

compute_viewshed <- function(dsm,
                             viewpoints,
                             offset_viewpoint=1.7,
                             offset_height = 0,
                             r = NULL,
                             refraction_factor = 0.13,
                             method = "plane",
                             parallel = FALSE,
                             workers = 1,
                             raster = FALSE,
                             plot = FALSE){
  multiviewpoints <- FALSE
  if (missing(dsm)) {
    stop("DSM is missing!")
  } else if (missing(viewpoints)) {
    stop("viewpoint(s) is missing!")
  }
  # Convert metre-based inputs to the DSM's native CRS unit (e.g. feet)
  if (is.null(r)) r <- 1000          # default radius in metres
  r                <- m_to_crs_units(r,                dsm)
  offset_viewpoint <- m_to_crs_units(offset_viewpoint, dsm)
  offset_height    <- m_to_crs_units(offset_height,    dsm)
  ext_dsm <- terra::ext(dsm)
  # if (dsm_units == "ft" && r > 3281) {
  #   r <- 3281
  # } else if (dsm_units == "m" && r > 1000) {
  #   r <- 1000
  # }
  if (plot) {
    raster <- TRUE
  }
  if (!(class(viewpoints)[1] == "numeric")) {
    if (class(viewpoints)[1] == "sf") {
      viewpoints <- sf::st_coordinates(viewpoints)
    } else {
      stop("If input viewpoints is not a vector or matrix, it has to be sf point(s)")
    }
  }
  if (is.vector(viewpoints)){
    viewpoints <- matrix(data = viewpoints,
                         nrow = 1,
                         ncol = 2)
  } else {
    if (length(viewpoints[,1]) > 1) {
      multiviewpoints <- TRUE
    }
  }
  if (multiviewpoints == FALSE){
    # validate the viewpoint
    if (viewpoints[,1] < ext_dsm[1] | viewpoints[,1] > ext_dsm[2]) {
      stop("The input viewpoint(s) is not on input dsm")
    }
    if (viewpoints[,2] < ext_dsm[3] | viewpoints[,2] > ext_dsm[4]) {
      stop("The input viewpoint(s) is not on input dsm")
    }
    viewpoints <- c(viewpoints[,1], viewpoints[,2])
    # compute viewshed
    output <- radius_viewshed(dsm, r,
                              refraction_factor,
                              viewpoints,
                              offset_viewpoint,
                              offset_height,
                              method)
    if (raster) {
      raster_data <- filter_invisible(output, raster)
      if (plot) {
        terra::plot(raster_data,
                    axes=FALSE,
                    box=FALSE,
                    legend = FALSE)
        v<- matrix(0,1,3)
        v[1,1] <- viewpoints[1]
        v[1,2] <- viewpoints[2]
        terra::plot(sp::SpatialPoints(v),
                    add=TRUE,
                    col="red",
                    axes=FALSE,
                    box=FALSE,
                    legend=FALSE)
      }
      return(raster_data)
    } else {
      return(output)
    }
  }else if (multiviewpoints){
    # validate the viewpoints
    if (min(viewpoints[,1]) < ext_dsm[1] | max(viewpoints[,1]) > ext_dsm[2]) {
      stop("Some of input viewpoint(s) is not on input dsm")
    }
    if (min(viewpoints[,2]) < ext_dsm[3] | max(viewpoints[,2]) > ext_dsm[4]) {
      stop("Some of input viewpoint(s) is not on input dsm")
    }
    inputs <- split(viewpoints,seq(nrow(viewpoints)))
    if (parallel == TRUE){
      if (workers > parallel::detectCores()) {
        workers <- parallel::detectCores()
      }
      # inputs <- split(viewpoints,seq(nrow(viewpoints)))
      if (isTRUE(Sys.info()[1]=="Windows") == TRUE){
        workers = 1
      }
      suppressWarnings(
        viewsheds <- paral_nix(X = inputs,
                               dsm = dsm,
                               r = r,
                               refraction_factor = refraction_factor,
                               offset = offset_viewpoint,
                               method = method,
                               workers = workers)
      )
    } else {
      suppressWarnings(
        viewsheds <- paral_nix(X = inputs,
                               dsm = dsm,
                               r = r,
                               refraction_factor = refraction_factor,
                               offset = offset_viewpoint,
                               method = method,
                               workers = 1)
      )
    }
    return(viewsheds)
  }
}
