test_that("fill_basins works", {
  filepath <- system.file("extdata", "terra_dem.tif", package = "flowdem")
  terra_dem <- terra::rast(filepath)
  u <- fill_basins(terra_dem)
  expect_equal(terra::as.matrix(terra::rast(filled), wide = TRUE),
             terra::as.matrix(u$dem, wide = TRUE))
})
