% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_diversity.R
\name{calculate_diversity}
\alias{calculate_diversity}
\title{calculate_diversity}
\usage{
calculate_diversity(viewshed, land, proportion = FALSE)
}
\arguments{
\item{viewshed}{Viewshed object.}

\item{land}{Raster. The raster of land use/land cover representing different
land use/cover classes.}

\item{proportion}{logical (Optional), indicating whether to return class
proportions along with the Shannon Diversity Index (SDI). (default is FALSE).}
}
\value{
List. a list containing the Shannon Diversity Index (SDI) and,
if the proportion parameter is set to TRUE, a table of class proportions
within the viewshed.
}
\description{
The calculate_diversity function is designed to calculate landscape
diversity metrics within a viewshed. It takes as input a land cover raster,
a viewshed object representing the observer's line of sight, and an optional
parameter to compute class proportions.
}
\examples{
library(viewscape)
# Load a viewpoint
test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
# load dsm raster
dsm <- terra::rast(system.file("test_dsm.tif", package ="viewscape"))
#Compute viewshed
output <- viewscape::compute_viewshed(dsm = dsm,
                                      viewpoints = test_viewpoint,
                                      offset_viewpoint = 6, r = 1600)
# load landuse raster
test_landuse <- terra::rast(system.file("test_landuse.tif",
                                        package ="viewscape"))
diversity <- viewscape::calculate_diversity(output,
                                            test_landuse)

}
