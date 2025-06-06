% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sector_mask.R
\name{fov_mask}
\alias{fov_mask}
\title{fov_mask}
\usage{
fov_mask(viewshed, fov)
}
\arguments{
\item{viewshed}{viewshed object, generated by compute_viewshed()}

\item{fov}{Vector, specifying the field of view with two angles in degree (e.g. c(10,100))
for masking a viewshed based on its viewpoints.
See details.}
}
\value{
viewshed object
}
\description{
The fov_mask function is designed to subset a viewshed based on
its viewpoint and the field of view
}
\details{
For defining the field of view ('fov'), angles range from 0 to 360 degrees,
with 0 inclusive and 360 exclusive. The initial angle must be smaller than the terminal
angle in the sequence c(a,b) (a < b). To capture the northeast quadrant of a viewshed,
one would use c(0,90), while the eastern quadrant would be delineated by c(45,315) as
shown below:
\tabular{rrr}{
135 \tab 90 \tab 45 \cr
180 \tab v \tab 0 \cr
225 \tab 270 \tab 315 \cr
}
Here, 'v' represents the viewpoint, with angles measured counterclockwise from due north.
}
\examples{
# Load a viewpoint
test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
# load dsm raster
dsm <- terra::rast(system.file("test_dsm.tif", package ="viewscape"))
# Compute viewshed
viewshed <- viewscape::compute_viewshed(dsm,
                                        viewpoints = test_viewpoint,
                                        offset_viewpoint = 6)
# subset viewshed using the field of view
out <- viewscape::fov_mask(viewshed, c(40,160))

}
\seealso{
\code{\link[=compute_viewshed]{compute_viewshed()}}
}
