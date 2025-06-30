<!-- badges: start -->
[![](https://www.r-pkg.org/badges/version/viewscape)](https://www.r-pkg.org/pkg/viewscape)
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![R-CMD-check](https://github.com/billbillbilly/viewscape/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/billbillbilly/viewscape/actions/workflows/R-CMD-check.yaml)
[![](https://cranlogs.r-pkg.org/badges/viewscape)](https://CRAN.R-project.org/package=viewscape)
![total](https://cranlogs.r-pkg.org/badges/grand-total/viewscape)
<!-- badges: end -->

# viewscape

<p align="left">

<img src="viewscape_hex-01.png" height="200">

</p>

## Introduction

The goal of `viewscape` package is to provide an accessible method of
carrying out landscape spatial analysis based on the viewshed within the
R environment. The viewscape R package can currently be installed via
github.

## Intallation

``` r

# Install the package from CRAN
install.packages("viewscape")

# Install the package from GitHub
devtools::install_github("land-info-lab/viewscape", dependencies=TRUE)

# Load the package
library(viewscape)
```

## Functionality and metrics

The basic viewshed analysis can be accessed by calling the
`compute_viewshed`. The two needed objects are a digital surface model
(DSM) and a viewpoint. It provides flexibility for single or
multi-viewpoint analyses and allows options for parallel processing,
raster output, and plotting. In addition, a `fov_mask` function is 
designed to subset a viewshed based on its viewpoint and the field of view. 

Based on the viewshed, visual magnitude, a gradient of visibility across the range, 
can be calculated using the `visual_magnitude` function. Also, a set of configuration metrics can be calculated using `calculate_viewmetrics`, `calculate_diversity`, and
`calculate_feature`.

The metrics are including:

-   Number of patches: Visible fragmentation measured by the total visible patches
    in the viewscape.
-   Mean shape index: Visible patchiness based on average perimeter-to-area ratio
    for all viewscape patches (vegetation and building). 
-   Edge density: A measure of visible complexity based on the length of
    patch edges per unit area.
-   Patch size: Total average size of a patch over the entire viewscape area.
-   Patch density: Visible landscape granularity based on measuring patch density.
-   Extent: The total area of the viewshed, calculated as the number of
    visible grid cells multiplied by the grid resolution.
-   Depth: The furthest visible distance within the viewshed from the
    viewpoint.
-   Vdepth: The standard deviation of distances to visible points,
    providing a measure of the variation in visible distances.
-   Horizontal: The total visible horizontal or terrestrial area within
    the viewshed.
-   Relief: The standard deviation of elevations of the visible ground
    surface.
-   Skyline: Variation of (Standard deviation) of the vertical viewscape
    (visible canopy and buildings).
-   Shannon diversity index: Based on the number of land use/cover
    classes and the proportion of distribution.
-   Proportion of other objects: Building, trees, or paved surface.



``` r
#Load in DSM
test_dsm <- terra::rast(system.file("test_dsm.tif", 
                                    package ="viewscape"))

#Load in the viewpoint
test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", 
                                          package = "viewscape"))

#Compute viewshed
output <- viewscape::compute_viewshed(dsm = test_dsm, 
                                      viewpoints = test_viewpoint, 
                                      offset_viewpoint = 6, 
                                      plot=TRUE)

# Load DTM
test_dtm <- terra::rast(system.file("test_dtm.tif", 
                                    package ="viewscape"))

# load landuse raster
test_landuse <- terra::rast(system.file("test_landuse.tif",
                                        package ="viewscape"))

# Load canopy raster
test_canopy <- terra::rast(system.file("test_canopy.tif", 
                                       package ="viewscape"))

# Load building footprints raster
test_building <- terra::rast(system.file("test_building.tif", 
                                       package ="viewscape"))

# calculate metrics given the viewshed
test_metrics <- viewscape::calculate_viewmetrics(output, 
                                                 test_dsm, 
                                                 test_dtm, 
                                                 list(test_canopy, test_building))

# the Shannon Diversity Index (SDI)
test_diversity <- calculate_diversity(test_landuse, output, proportion = TRUE)
```

From viewshed analysis, the visible area of a viewpoint is presented by
visible points. Several viewshed metrics can be
calculated based on the visible points. For further information on these
metrics and the rest of the functions available in this package, please
refer to the [package
website](https://land-info-lab.github.io/viewscape/). For more
information and examples of the functions, check out the [package
vignette](https://github.com/land-info-lab/viewscape/blob/master/viewscape.md).

## Note
This package may take a long time to run if using spatially large or
high-resolution digital elevation models.

The package currently does not support multi-core processing on Windows system. 
The function 'compute_viewshed' with 'parallel = TRUE' will automatically work
on single worker.

## Reference

Franklin, W. R., & Ray, C. (1994, May).
Higher isn’t necessarily better: Visibility algorithms and experiments.
In Advances in GIS research: sixth international symposium on spatial
data handling (Vol. 2, pp. 751-770). Edinburgh: Taylor & Francis.

Wang, J., Robinson, G. J., & White, K. (2000).
Generating viewsheds without using sightlines.
Photogrammetric engineering and remote sensing, 66(1), 87-90.

Chamberlain, B. C., & Meitner, M. J. (2013). 
A route-based visibility analysis for landscape management. 
Landscape and Urban Planning, 111, 13-24.

Tabrizian, P., Baran, P. K., Van Berkel, D., Mitasova, H., & Meentemeyer, R. (2020). Modeling restorative potential of urban environments by coupling viewscape analysis of lidar data with experiments in immersive virtual environments. Landscape and Urban planning, 195, 103704.

## Issues and bugs

If you discover a bug that is not already a [reported
issue](https://github.com/billbillbilly/viewscape/issues), please [open
a new issue](https://github.com/billbillbilly/viewscape/issues/new)
providing a reproducible example.
