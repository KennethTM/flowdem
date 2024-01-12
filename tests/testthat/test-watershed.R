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
  
  # Load poly
  filepath <- system.file("extdata", "poly.gpkg", package = "flowdem")
  poly <- sf::st_read(filepath, quiet = TRUE)
  
  # Load watershed of poly
  filepath <- system.file("extdata", "watershed_poly.tif", package = "flowdem")
  expected_watershed_poly <- terra::rast(filepath)
  
  # Load poly raster
  filepath <- system.file("extdata", "poly_rast.tif", package = "flowdem")
  poly_rast <- terra::rast(filepath)
  
  # Load watershed of poly raster
  filepath <- system.file("extdata", "watershed_poly_rast.tif", package = "flowdem")
  expected_watershed_poly_rast <- terra::rast(filepath)

  # Test watershed of point and of line
  actual_watershed_point <- watershed(dirs, p)
  expect_equal(terra::compareGeom(expected_watershed_point, actual_watershed_point), TRUE)
  expect_equal(unname(terra::values(expected_watershed_point)), unname(terra::values(actual_watershed_point)))

  actual_watershed_line <- watershed(dirs, l)
  expect_equal(terra::compareGeom(expected_watershed_line, actual_watershed_line), TRUE)
  expect_equal(unname(terra::values(expected_watershed_line)), unname(terra::values(actual_watershed_line)))

  # Test watershed of polygon
  actual_watershed_poly <- watershed(dirs, poly)
  expect_equal(terra::compareGeom(expected_watershed_poly, actual_watershed_poly), TRUE)
  expect_equal(unname(terra::values(expected_watershed_poly)), unname(terra::values(actual_watershed_poly)))
  
  # Test watershed of raster of polygon
  actual_watershed_poly_rast <- watershed(dirs, poly_rast)
  expect_equal(terra::compareGeom(expected_watershed_poly_rast, actual_watershed_poly_rast), TRUE)
  expect_equal(unname(terra::values(expected_watershed_poly_rast)), unname(terra::values(actual_watershed_poly_rast)))
  
  # Test that watershed computed using polygon and raster of same polygon is identical
  expect_equal(terra::compareGeom(actual_watershed_poly, actual_watershed_poly_rast), TRUE)
  expect_equal(unname(terra::values(actual_watershed_poly)), unname(terra::values(actual_watershed_poly_rast)))
  
})
