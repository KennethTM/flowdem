test_that("fill_basins works", {

  # Load original DEM
  filepath <- system.file("extdata", "dem.tif", package = "flowdem")
  dem <- terra::rast(filepath)

  # Load filled DEM
  filepath <- system.file("extdata", "filled.tif", package = "flowdem")
  expected_filled <- terra::rast(filepath)

  # Load basins
  filepath <- system.file("extdata", "basins.tif", package = "flowdem")
  expected_basins <- terra::rast(filepath)

  # Test fill
  actual <- fill_basins(dem)
  actual_filled <- actual$dem
  actual_basins <- actual$basins

  expect_equal(terra::compareGeom(expected_filled, actual_filled), TRUE)
  expect_equal(unname(terra::values(expected_filled)), unname(terra::values(actual_filled)))

  expect_equal(terra::compareGeom(expected_basins, actual_basins), TRUE)
  expect_equal(unname(terra::values(expected_basins)), unname(terra::values(actual_basins)))
})
