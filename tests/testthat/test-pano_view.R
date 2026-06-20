testthat::test_that("pano_view returns a panorama raster", {
  test_dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
  test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp",
                                            package = "viewscape"))

  result <- viewscape::pano_view(test_dsm, test_viewpoint, 6,
                                 method = "cylindrical",
                                 pano_dim = c(12, 24),
                                 max_dist = 20,
                                 step_size = 2)

  testthat::expect_s4_class(result, "SpatRaster")
  testthat::expect_equal(dim(result), c(12, 24, 1))
  testthat::expect_true(any(!is.na(terra::values(result))))
})

testthat::test_that("pano_view supports equirectangular method", {
  test_dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
  test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp",
                                            package = "viewscape"))

  result <- viewscape::pano_view(dsm = test_dsm,
                                 vpt = test_viewpoint,
                                 method = "equirectangular",
                                 pano_dim = c(8, 16),
                                 max_dist = 20,
                                 step_size = 2)

  testthat::expect_s4_class(result, "SpatRaster")
  testthat::expect_equal(dim(result), c(8, 16, 1))
})

testthat::test_that("pano_view uses distinct cylindrical and equirectangular projections", {
  test_dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
  test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp",
                                            package = "viewscape"))

  cylindrical <- viewscape::pano_view(dsm = test_dsm,
                                      vpt = test_viewpoint,
                                      method = "cylindrical",
                                      pano_dim = c(12, 24),
                                      max_dist = 30,
                                      step_size = 2)
  equirectangular <- viewscape::pano_view(dsm = test_dsm,
                                          vpt = test_viewpoint,
                                          method = "equirectangular",
                                          pano_dim = c(12, 24),
                                          max_dist = 30,
                                          step_size = 2)

  testthat::expect_false(isTRUE(all.equal(terra::values(cylindrical),
                                          terra::values(equirectangular),
                                          check.attributes = FALSE)))
})

testthat::test_that("pano_view returns semantic panorama when semantic is supplied", {
  test_dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
  test_semantic <- terra::rast(system.file("test_landuse.tif", package = "viewscape"))
  test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp",
                                            package = "viewscape"))

  testthat::expect_warning(
    result <- viewscape::pano_view(dsm = test_dsm,
                                   vpt = test_viewpoint,
                                   semantic = test_semantic,
                                   method = "equirectangular",
                                   pano_dim = c(8, 16),
                                   max_dist = 20,
                                   step_size = 2),
    "resampling semantic to the dsm grid"
  )

  testthat::expect_s4_class(result, "SpatRaster")
  testthat::expect_equal(dim(result), c(8, 16, 2))
  testthat::expect_equal(names(result),
                         c("equirectangular_depth", "equirectangular_semantic"))
  testthat::expect_true(any(!is.na(terra::values(result[[2]]))))
})

testthat::test_that("equirectangular semantic projection extrudes binary masks", {
  test_dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
  test_building <- terra::rast(system.file("test_building.tif",
                                           package = "viewscape"))
  test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp",
                                            package = "viewscape"))

  testthat::expect_warning(
    result <- viewscape::pano_view(dsm = test_dsm,
                                   vpt = test_viewpoint,
                                   semantic = test_building,
                                   method = "equirectangular",
                                   pano_dim = c(32, 64),
                                   max_dist = 300,
                                   step_size = 4),
    "resampling semantic to the dsm grid"
  )

  semantic_matrix <- terra::as.matrix(result[[2]], wide = TRUE)
  occupied_rows <- which(rowSums(semantic_matrix > 0.5, na.rm = TRUE) > 0)

  testthat::expect_gt(sum(semantic_matrix > 0.5, na.rm = TRUE), 500)
  testthat::expect_gt(length(occupied_rows), 10)
  testthat::expect_equal(range(which(colSums(semantic_matrix > 0.5,
                                             na.rm = TRUE) > 0)),
                         c(1, 64))
})

testthat::test_that("pano_view validates required inputs", {
  test_dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))

  testthat::expect_error(viewscape::pano_view(dsm = test_dsm), "vpt is missing")
  testthat::expect_error(viewscape::pano_view(dsm = test_dsm,
                                              vpt = c(1)),
                         "vpt must contain x and y coordinates")
  testthat::expect_error(viewscape::pano_view(dsm = test_dsm,
                                              vpt = c(1, 1),
                                              step_size = 0),
                         "step_size")
  testthat::expect_error(viewscape::pano_view(dsm = test_dsm,
                                              vpt = c(1, 1),
                                              semantic = "landuse"),
                         "semantic must be a SpatRaster")
  bad_semantic <- terra::rast(nrows = 2, ncols = 2)
  testthat::expect_error(viewscape::pano_view(dsm = test_dsm,
                                              vpt = c(1, 1),
                                              semantic = bad_semantic),
                         "semantic has no values")
})
