test_that("accum works", {

  # Load d8 dirs
  filepath <- system.file("extdata", "dirs.tif", package = "flowdem")
  dirs <- terra::rast(filepath)

  # Load accum
  filepath <- system.file("extdata", "accum.tif", package = "flowdem")
  expected_accum <- terra::rast(filepath)

  # Test fill
  actual_accum <- accum(dirs)
  expect_equal(terra::compareGeom(expected_accum, actual_accum), TRUE)
  expect_equal(unname(terra::values(expected_accum)), unname(terra::values(actual_accum)))
})
