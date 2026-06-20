testthat::test_that("returns matrix of correct dimensions", {

  test_dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
  test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))

  # Duplicate the viewpoint to create a 2-point set for testing
  vp2 <- rbind(sf::st_coordinates(test_viewpoint),
               sf::st_coordinates(test_viewpoint) + 10)

  mat <- viewscape::intervis_network(vp2, test_dsm, offset_viewpoint = 1.7)

  testthat::expect_true(is.matrix(mat))
  testthat::expect_equal(dim(mat), c(2L, 2L))
  testthat::expect_true(is.na(mat[1, 1]))
  testthat::expect_true(is.na(mat[2, 2]))
  testthat::expect_true(mat[1, 2] %in% c(0L, 1L))
  testthat::expect_true(mat[2, 1] %in% c(0L, 1L))
})

testthat::test_that("lines output returns sf object", {

  test_dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
  test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))

  vp2 <- rbind(sf::st_coordinates(test_viewpoint),
               sf::st_coordinates(test_viewpoint) + 10)

  net <- viewscape::intervis_network(vp2, test_dsm, offset_viewpoint = 1.7,
                                     output = "lines")

  testthat::expect_s3_class(net, "sf")
  testthat::expect_true(all(c("from", "to", "mutual") %in% names(net)))
})

testthat::test_that("per-viewpoint heights accepted", {

  test_dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
  test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))

  vp2 <- rbind(sf::st_coordinates(test_viewpoint),
               sf::st_coordinates(test_viewpoint) + 10)

  testthat::expect_no_error(
    viewscape::intervis_network(vp2, test_dsm, offset_viewpoint = c(1.7, 6))
  )
})

testthat::test_that("invalid inputs are caught", {

  test_dsm <- terra::rast(system.file("test_dsm.tif", package = "viewscape"))
  test_viewpoint <- sf::read_sf(system.file("test_viewpoint.shp", package = "viewscape"))
  vp2 <- rbind(sf::st_coordinates(test_viewpoint),
               sf::st_coordinates(test_viewpoint) + 10)

  testthat::expect_error(viewscape::intervis_network(vp2, test_dsm, output = "bad"),
                         "output must be")
  testthat::expect_error(viewscape::intervis_network(vp2, test_dsm,
                                                     offset_viewpoint = c(1, 2, 3)),
                         "offset_viewpoint")
  testthat::expect_error(viewscape::intervis_network(missing, test_dsm),
                         "viewpoints is missing")
})
