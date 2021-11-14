
<!-- badges: start -->

[![R-CMD-check](https://github.com/billbillbilly/viewscape/workflows/R-CMD-check/badge.svg)](https://github.com/billbillbilly/viewscape/actions)
[![Codecov test
coverage](https://codecov.io/github/billbillbilly/viewscape/branch/master/graph/badge.svg)](https://codecov.io/github/billbillbilly/viewscape?branch=master)
<!-- badges: end -->

# viewscape

<p align="left">

<img src=".//man//figures//viewscape_hex.png" height="200">

</p>

## Introduction

The goal of viewscape package is to provide an accessible method of
carrying out viewshed analysis within the R environment. The viewscape R
pacakge can currently be downloaded via github.

``` r
library(devtools)

install_github("billbillbilly/viewscape")
```

The basic viewshed analysis can be accessed through calling the
`calculate_viewshed`. The two needed objects to calculate the viewshed
are a digital surface model (DSM) and a viewpoint.

``` r
  #Load in DSM
  test_dsm <- raster::raster(system.file("test_data\\test_dsm.tif",
                                         package ="viewscape"))

  #Load in the viewpoint
  test_viewpoint <- sf::read_sf(system.file("test_data\\test_viewpoint.shp",
                                               package = "viewscape"))

  #Transform viewpoint from shape file to coordinates 
  test_viewpoint <- sf::st_coordinates(test_viewpoint)
  test_viewpoint <- c(test_viewpoint[,1], test_viewpoint[,2])

  #Run function
  test_visiblepoint <- viewscape::calculate_viewshed(dsm = test_dsm,
                                                     viewpoint = test_viewpoint)
```

From viewshed analysis, the visible area of a viewpoint is presented by visible points. There are several viewshed metrics such as can be calculated based on the visible points.

The function of view depth analysis can calculate two different metrics: the furthest distance and standard deviation of distances. To calculate view depth, there are two needed objects: the DSM that was used to get viewshed and result from viewshed analysis. Additionally, the third input parameter is used to customize the output.

``` r
   # calculate view depth within the viewshed
  test_depth <- viewscape::get_depth(test_visiblepoint, test_viewpoint, 3)
```

The function of extent analysis can calculate the total area of viewshed and needs the DSM that was used to get viewshed and result from viewshed analysis. 

``` r
  # calculate extent of the viewshed
  test_extent <- viewscape::get_extent(test_dsm, test_visiblepoint)
```

The following function can calculate the area of ground surface and standard deviation of elevations within a viewshed. The function needs a DSM and a DEM/DTM to calculate the metrics. Additionally, the third input parameter is used to customize the output.

``` r
  # load DEM
  test_dem <- raster::raster(system.file("test_data/test_dem.tif",
                                         package ="viewscape"))
  # calculate the area of ground surface and standard deviation of elevations
   test_horizontal <- viewscape::get_horizontal(test_dsm, test_dem, 
                                                test_visiblepoint, 3)
```

To calculate canopy area in a viewshed, the DSM that was used to get viewshed and a raster of canopy are needed. Additionally, input parameter 'data' is used to indicate the type of input canopy raster and 'nodata' is used to indicate the the value of cells that don't have any canopy.

``` r
  # load canopy raster
  test_canopy <- raster::raster(system.file("test_data/test_canopy.tif",
                                            package ="viewscape"))
  # calculate the area of canopy
  test_canopy_area <- viewscape::calculate_canopy(data = 1, test_canopy, 
                                                  nodata = 0, test_dsm, 
                                                  test_visiblepoint)
```

The land cover analysis calculates the areas of perviousness and imperviousness and percentages(%) of perviousness and imperviousness in a viewshed. The DSM that was used to get viewshed, visible points, and a raster of land cover are needed. Input parameters 'vegetation' and 'imperviousness' are used to indicate the code of perviousness including trees or grass and the code of imperviousness including buildings, parking, and roads.

``` r
  # load landuse raster
 test_landcover <- raster::raster(system.file("test_data/test_landcover.tif",
                                            package ="viewscape"))
  # calculate the areas and percentages of perviousness and imperviousness 
  # in the sample data of land cover, value 2 is for vegetation and value 4 is for imperviousness
  test_landcover_area <- viewscape::calculate_landcover(landcover = test_landcover, 
                                                        dsm = test_dsm,
                                                        visiblepoints = test_visiblepoint,
                                                        vegetation = 2, imperviousness = 4)
```

The land use analysis calculates the percentage(%) of each type of land use in a viewshed. The DSM that was used to get viewshed, visible points, and a raster of land use are needed.

``` r
  # load landuse raster
 test_landuse <- raster::raster(system.file("test_data/test_landuse.tif",
                                            package ="viewscape"))
  # calculate the percentage of each type of land use  
  test_landuse_area <- viewscape::calculate_landuse(landuse = test_landuse, 
                                                    dsm = test_dsm,
                                                    visiblepoints = test_visiblepoint)
```

For further information on the rest of the functions available in this
package please refer to the [package
website](https://billbillbilly.github.io/viewscape/). For more
information and examples of the functions check out the [package
vignette](needs%20to%20be%20created).

## Issues and bugs

This package may take a long time to run if using spatially large or
high resolution digital elevation models.

If you discover a bug not associated with connection to the API that is
not already a [reported
issue](https://github.com/billbillbilly/viewscape/issues), please [open
a new issue](https://github.com/billbillbilly/viewscape/issues/new)
providing a reproducible example.
