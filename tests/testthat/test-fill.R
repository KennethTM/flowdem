test_that("fill works", {
  expect_equal(raster::values(filled),
               raster::values(fill(test_dem, 0.01)))
})
