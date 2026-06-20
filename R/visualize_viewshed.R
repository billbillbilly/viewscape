#' visualize_viewshed
#' @description The visualize_viewshed function is designed for the visualization
#' of a viewshed analysis, providing users with various options for visualizing
#' the results. The function works with a viewshed object and offers multiple
#' plotting and output types.
#'
#' @param viewshed Viewshed object
#' @param plottype Character, specifying the type of visualization ("polygon" or
#' "raster").
#' @param outputtype Character, specifying the type of output object ("raster"
#' or "polygon").
#' @return Visualized viewshed as either a raster or polygon object,
#' depending on the outputtype specified.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Load a viewpoint
#' test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
#' # load dsm raster
#' dsm <- terra::rast(system.file("test_dsm.tif", package ="viewscape"))
#' #Compute viewshed
#' viewshed <- compute_viewshed(dsm = dsm,
#'                              viewpoints = test_viewpoint,
#'                              offset_viewpoint = 6)
#' # Visualize the viewshed as polygons
#' visualize_viewshed(viewshed, plottype = "polygon")
#' # Visualize the viewshed as a raster
#' visualize_viewshed(viewshed, plottype = "raster")
#' # Get the visualized viewshed as a polygon object
#' polygon_viewshed <- visualize_viewshed(viewshed,
#'                                        plottype = "polygon",
#'                                        outputtype = "polygon")
#'}

visualize_viewshed <- function(viewshed,
                               plottype = "",
                               outputtype = "") {
  if (missing(viewshed)){
    stop("Viewshed object is missing")
  }
  valid_plottypes   <- c("", "polygon", "raster")
  valid_outputtypes <- c("", "raster", "polygon")
  if (!plottype %in% valid_plottypes) {
    warning("Unrecognized plottype '", plottype, "'. Use 'polygon' or 'raster'.")
  }
  if (!outputtype %in% valid_outputtypes) {
    warning("Unrecognized outputtype '", outputtype, "'. Use 'raster' or 'polygon'.")
  }

  # vectorize the viewshed
  mask_v <- get_patch(viewshed)

  # pre-compute polygon once if needed by either plottype or outputtype
  if (plottype == "polygon" || outputtype == "polygon") {
    polygon_v <- terra::as.polygons(mask_v)
  }

  if (plottype == "polygon"){
    terra::plot(polygon_v, col = rgb(0, 1, 0, 0.3), border = NA)
  } else if (plottype == "raster"){
    terra::plot(mask_v)
  }

  if (outputtype == "raster"){
    return(mask_v)
  } else if (outputtype == "polygon"){
    return(sf::st_as_sf(polygon_v))
  }

  invisible(NULL)
}
