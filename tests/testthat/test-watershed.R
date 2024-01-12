test_that("watershed works", {

  # Load d8 dirs
  filepath <- system.file("extdata", "dirs.tif", package = "flowdem")
  dirs <- terra::rast(filepath)

  # Load point
  filepath <- system.file("extdata", "point.gpkg", package = "flowdem")
  p <- sf::st_read(filepath, quiet = TRUE)

  # Load watershed of point
  filepath <- system.file("extdata", "watershed_point.tif", package = "flowdem")
  expected_watershed_point <- terra::rast(filepath)

  # Load line
  filepath <- system.file("extdata", "line.gpkg", package = "flowdem")
  l <- sf::st_read(filepath, quiet = TRUE)

  # Load watershed of line
  filepath <- system.file("extdata", "watershed_line.tif", package = "flowdem")
  expected_watershed_line <- terra::rast(filepath)

  # Test watershed of point and of line
  actual_watershed_point <- watershed(dirs, p)
  expect_equal(terra::compareGeom(expected_watershed_point, actual_watershed_point), TRUE)
  expect_equal(unname(terra::values(expected_watershed_point)), unname(terra::values(actual_watershed_point)))

  actual_watershed_line <- watershed(dirs, l)
  expect_equal(terra::compareGeom(expected_watershed_line, actual_watershed_line), TRUE)
  expect_equal(unname(terra::values(expected_watershed_line)), unname(terra::values(actual_watershed_line)))

})
