test_that("breach works", {

  # Load original DEM
  filepath <- system.file("extdata", "dem.tif", package = "flowdem")
  dem <- terra::rast(filepath)

  # Load breached DEM
  filepath <- system.file("extdata", "breached.tif", package = "flowdem")
  expected_breached <- terra::rast(filepath)

  # Test fill
  actual_breached <- breach(dem)
  expect_equal(terra::compareGeom(expected_breached, actual_breached), TRUE)
  expect_equal(unname(terra::values(expected_breached)), unname(terra::values(actual_breached)))

})
