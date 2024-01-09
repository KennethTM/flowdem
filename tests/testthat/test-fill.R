test_that("fill works", {
  filepath <- system.file("extdata", "terra_dem.tif", package = "flowdem")
  terra_dem <- terra::rast(filepath)
    expect_equal(terra::as.matrix(terra::rast(filled), wide = TRUE),
               terra::as.matrix(fill(terra_dem, 0.01), wide = TRUE))
})