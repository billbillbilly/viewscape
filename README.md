<!-- badges: start -->
[![](https://www.r-pkg.org/badges/version/viewscape)](https://www.r-pkg.org/pkg/viewscape)
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![R-CMD-check](https://github.com/land-info-lab/viewscape/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/land-info-lab/viewscape/actions/workflows/R-CMD-check.yaml)
[![](https://cranlogs.r-pkg.org/badges/viewscape)](https://CRAN.R-project.org/package=viewscape)
![total](https://cranlogs.r-pkg.org/badges/grand-total/viewscape)
<!-- badges: end -->

# viewscape

<p align="left">

<img src="viewscape_hex-01.png" height="200">

</p>

## Introduction

The goal of viewscape is to provide accessible landscape spatial analysis
based on the viewshed within the R environment. The package computes
viewsheds from digital surface models, derives visual magnitude and a
rich set of viewscape configuration metrics, generates panoramic views,
and maps intervisibility networks between viewpoints.

## Installation

``` r
# Install the package from CRAN
install.packages("viewscape")

# Install the development version from GitHub
devtools::install_github("land-info-lab/viewscape", dependencies = TRUE)

# Load the package
library(viewscape)
```

## Functions

| Function | Description |
|---|---|
| `compute_viewshed()` | Compute a viewshed from a DSM and one or more viewpoints. Supports single and multi-viewpoint analysis, parallel processing, and raster output. Observer height (`offset_viewpoint`), target height (`offset_height`), and radius (`r`) are all specified **in metres** and converted automatically to the DSM’s CRS unit. |
| `visual_magnitude()` | Calculate visual magnitude — the gradient of visibility across the viewshed — following Chamberlain & Meitner (2013). |
| `calculate_viewmetrics()` | Compute a set of configuration metrics (extent, depth, relief, skyline, patch statistics) for a viewshed. |
| `calculate_diversity()` | Compute the Shannon Diversity Index of land cover within a viewshed. |
| `calculate_feature()` | Compute the proportion of a specific land-use or land-cover class visible in a viewshed. |
| `sector_mask()` | Subset a viewshed to a field-of-view sector defined by a bearing and angular width. |
| `intervis_network()` | Compute pairwise intervisibility among a set of viewpoints and return results as an adjacency matrix or an sf network of lines. |
| `pano_view()` | Generate a cylindrical or equirectangular panoramic view from a DSM. Supports satellite RGB colour (via greenSD), semantic land-cover layers, and customisable ray-marching parameters. |
| `visualize_viewshed()` | Visualize a viewshed object with various plotting options. |

## Quick start

``` r
# Load data
test_dsm      <- terra::rast(system.file("test_dsm.tif",      package = "viewscape"))
test_dtm      <- terra::rast(system.file("test_dtm.tif",      package = "viewscape"))
test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
test_canopy   <- terra::rast(system.file("test_canopy.tif",   package = "viewscape"))
test_building <- terra::rast(system.file("test_building.tif", package = "viewscape"))
test_landuse  <- terra::rast(system.file("test_landuse.tif",  package = "viewscape"))

# Compute viewshed (observer 1.7 m above surface, 1000 m radius — both in metres)
vs <- viewscape::compute_viewshed(dsm = test_dsm,
                                  viewpoints = test_viewpoint,
                                  offset_viewpoint = 1.7,
                                  r = 1000)

# Visual magnitude
vm <- viewscape::visual_magnitude(vs, test_dsm)

# Configuration metrics
metrics <- viewscape::calculate_viewmetrics(vs, test_dsm, test_dtm,
                                            list(test_canopy, test_building))

# Shannon Diversity Index
diversity <- viewscape::calculate_diversity(test_landuse, vs, proportion = TRUE)

# Panoramic view (cylindrical, no internet needed)
pano <- viewscape::pano_view(test_dsm, test_viewpoint, h = 1.7,
                             method = "cylindrical", plot = TRUE)

# Intervisibility network (returns adjacency matrix by default)
# mat <- viewscape::intervis_network(test_viewpoint, test_dsm)
```

## Viewscape metrics

`calculate_viewmetrics()` returns the following metrics:

- **Number of patches (Nump)**: visible fragmentation — total visible patches within the viewscape.
- **Mean shape index (MSI)**: visible patchiness based on average perimeter-to-area ratio for all viewscape patches.
- **Edge density (ED)**: visible complexity based on the length of patch edges per unit area.
- **Patch size (PS)**: total average size of a patch over the entire viewscape area.
- **Patch density (PD)**: visible landscape granularity.
- **Extent**: total visible area (visible cells × resolution²).
- **Depth**: furthest visible distance from the viewpoint.
- **Vdepth**: standard deviation of distances to visible points.
- **Horizontal**: total visible terrestrial (ground-level) area.
- **Relief**: standard deviation of visible ground surface elevations.
- **SVF (Sky View Factor)**: proportion of the visible sky hemisphere at the viewpoint, computed by casting rays in 36 azimuth directions and calculating SVF = mean(cos²(max obstruction angle)). Ranges from 0 (fully enclosed) to 1 (open sky). Method follows Oke (1981) as implemented in the shadow package (Dorman et al. 2019).
- **Skyline**: standard deviation of DSM heights for cells with visible canopy or buildings.

## Note

The package does not support multi-core processing on Windows.
`compute_viewshed(parallel = TRUE)` will automatically fall back to a
single worker on Windows.

## References

Franklin, W. R., & Ray, C. (1994). Higher isn’t necessarily better: Visibility algorithms and experiments. *Advances in GIS Research*, 2, 751–770.

Wang, J., Robinson, G. J., & White, K. (2000). Generating viewsheds without using sightlines. *Photogrammetric Engineering and Remote Sensing*, 66(1), 87–90.

Chamberlain, B. C., & Meitner, M. J. (2013). A route-based visibility analysis for landscape management. *Landscape and Urban Planning*, 111, 13–24.

Tabrizian, P., Baran, P. K., Van Berkel, D., Mitasova, H., & Meentemeyer, R. (2020). Modeling restorative potential of urban environments by coupling viewscape analysis of lidar data with experiments in immersive virtual environments. *Landscape and Urban Planning*, 195, 103704.

Oke, T. R. (1981). Canyon geometry and the nocturnal urban heat island: comparison of scale model and field observations. *Journal of Climatology*, 1(3), 237–254.

Dorman, M., Vulkan, A., Erell, E., & Kloog, I. (2019). shadow: Geometric Shadow Calculations. *The R Journal*, 11(1), 287–309. https://github.com/michaeldorman/shadow

## Issues and bugs

This package may take a long time to run with spatially large or
high-resolution elevation models.

If you discover a bug not already in the [reported
issues](https://github.com/land-info-lab/viewscape/issues), please [open
a new issue](https://github.com/land-info-lab/viewscape/issues/new)
with a reproducible example.
