% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_viewshed.R
\name{compute_viewshed}
\alias{compute_viewshed}
\title{compute_viewshed}
\usage{
compute_viewshed(
  dsm,
  viewpoints,
  offset_viewpoint = 1.7,
  offset_height = 0,
  r = NULL,
  refraction_factor = 0.13,
  method = "plane",
  parallel = FALSE,
  workers = 1,
  raster = FALSE,
  plot = FALSE
)
}
\arguments{
\item{dsm}{Raster, the digital surface model/digital elevation model}

\item{viewpoints}{sf point(s) or vector including x,y coordinates of a viewpoint
or a matrix including several viewpoints with x,y coordinates}

\item{offset_viewpoint}{numeric, setting the height of the viewpoint.
(default is 1.7 meters).}

\item{offset_height}{numeric, setting the height of positions that a given
viewpoint will look at. (defaut is 0)}

\item{r}{Numeric (optional), setting the radius for viewshed analysis.
(The default is 1000m/3281ft)}

\item{refraction_factor}{Number, indicating the refraction factor.
The refraction factor adjusts the effect of atmospheric refraction
on the apparent curvature of the Earth. In most standard applications, a refraction factor
of 0.13 is used, and so does this function. However, the appropriate refraction factor may
vary depending on environmental conditions.}

\item{method}{Character, The algorithm for computing a viewshed:
"plane" and "los" (see details). "plane" is used as default.}

\item{parallel}{Logical, (default is FALSE) indicating if parallel computing
should be used to compute viewsheds of multiview points. When it is TRUE,
arguements 'raster' and 'plot' are ignored}

\item{workers}{Numeric, indicating the number of CPU cores that will be used
for parallel computing. It is required if 'parallel' is 'TRUE'.}

\item{raster}{Logical, (default is FALSE) if it is TRUE, the raster of
viewshed will be returned.
The default is FALSE}

\item{plot}{Logical, (default is FALSE) if it is TRUE, the raster of
viewshed will be displayed}
}
\value{
Raster or list. For single-viewpoint analysis, the function returns
either a raster (raster is TRUE) or a viewshed object. Value 1 means visible while
value 0 means invisible. For multi-viewpoint analysis, a list of viewsheds
is returned.
}
\description{
The compute_viewshed function is designed for computing viewsheds,
which are areas visible from specific viewpoints, based on a Digital Surface
Model (DSM). It provides flexibility for single or multi-viewpoint analyses
and allows options for parallel processing, raster output, and plotting.
}
\details{
For method, "plane" is the reference plane algorithm introduced by
Wang et al. (2000) and "los" is the line of sight algorithm (Franklin & Ray, 1994).
The reference plane algorithm can be more time-efficient than the line of sight
algorithm, whereas the accuracy of the line of sight is better.
}
\examples{
# Load a viewpoint
test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
# load dsm raster
dsm <- terra::rast(system.file("test_dsm.tif", package ="viewscape"))
#Compute viewshed
output <- viewscape::compute_viewshed(dsm = dsm,
                                      viewpoints = test_viewpoint,
                                      offset_viewpoint = 6, r = 1600)
}
\references{
Franklin, W. R., & Ray, C. (1994, May).
Higher isn’t necessarily better: Visibility algorithms and experiments.
In Advances in GIS research: sixth international symposium on spatial
data handling (Vol. 2, pp. 751-770). Edinburgh: Taylor & Francis.

Wang, J., Robinson, G. J., & White, K. (2000).
Generating viewsheds without using sightlines.
Photogrammetric engineering and remote sensing, 66(1), 87-90.
}
\seealso{
\code{\link[=fov_mask]{fov_mask()}} \code{\link[=visual_magnitude]{visual_magnitude()}}
}
