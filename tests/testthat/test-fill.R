test_that("fill works", {
  # expect_equal(raster::values(filled),
  #              raster::values(fill(test_dem, 0.01)))
  filepath <- system.file("extdata", "terra_dem.tif", package = "flowdem")
  terra_dem <- terra::rast(filepath)
    expect_equal(raster::values(filled),
               terra::values(fill(terra_dem, 0.01)))
})
