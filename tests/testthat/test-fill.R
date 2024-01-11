test_that("fill works", {

  # Load original DEM
  filepath <- system.file("extdata", "dem.tif", package = "flowdem")
  dem <- terra::rast(filepath)

  # Load breached DEM
  filepath <- system.file("extdata", "breached.tif", package = "flowdem")
  breached <- terra::rast(filepath)

  # Load filled DEM without epsilon. This was generated from the original DEM.
  filepath <- system.file("extdata", "filled.tif", package = "flowdem")
  expected_filled <- terra::rast(filepath)

  # Test fill
  actual_filled <- fill(dem)
  expect_equal(terra::compareGeom(expected_filled, actual_filled), TRUE)
  expect_equal(unname(terra::values(expected_filled)), unname(terra::values(actual_filled)))

  # Load filled with epsilon DEM. This was generated from the breached DEM.
  filepath <- system.file("extdata", "filled_eps.tif", package = "flowdem")
  expected_filled <- terra::rast(filepath)

  # Test fill with epsilon
  actual_filled <- fill(breached, epsilon = 0.01)
  expect_equal(terra::compareGeom(expected_filled, actual_filled), TRUE)
  expect_equal(unname(terra::values(expected_filled)), unname(terra::values(actual_filled)))

})