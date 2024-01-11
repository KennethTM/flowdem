test_that("dirs works", {

  # Load breached and filled with epsilon DEM
  filepath <- system.file("extdata", "filled_eps.tif", package = "flowdem")
  filled_eps <- terra::rast(filepath)

  # Load d8 dirs
  filepath <- system.file("extdata", "dirs.tif", package = "flowdem")
  expected_dirs <- terra::rast(filepath)

  # Test dirs
  actual_dirs <- dirs(filled_eps)
  expect_equal(terra::compareGeom(actual_dirs, expected_dirs), TRUE)
  expect_equal(unname(terra::values(actual_dirs)), unname(terra::values(expected_dirs)))

})
