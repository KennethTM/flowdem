# This script documents the process of generating the extdata for the package

# dem.tif is extracted from the SRTM 1-arcsecond global DEM

# devpath is the path to the root of the flowdem package on your machine
devpath <- "/media/kenneth/d6c13395-8492-49ee-9c0f-6a165e34c95c1/flowdem"

# Read the DEM
dem <- terra::rast(paste0(devpath, "/inst/extdata/dem.tif"))

# Fill the DEM, leaving flat areas flat, and compute coastal basins simultaneously
filled_basins <- flowdem::fill_basins(dem)
terra::writeRaster(filled_basins$dem, paste0(devpath, "/inst/extdata/filled.tif"), overwrite = TRUE)
terra::writeRaster(filled_basins$basins, paste0(devpath, "/inst/extdata/basins.tif"), overwrite = TRUE)

# Breach the DEM, carving through obstacles to ensure flow
breached <- flowdem::breach(dem)
terra::writeRaster(breached, paste0(devpath, "/inst/extdata/breached.tif"), overwrite = TRUE)

# Fill the depressions in breached DEM, adding a small gradient on flat surfaces to ensure flow
# We set the datatype to FLT8S, see documentation of fill() to understand why
filled_eps <- flowdem::fill(breached, epsilon = 0.01)
terra::writeRaster(filled_eps, paste0(devpath, "/inst/extdata/filled_eps.tif"), overwrite = TRUE, datatype="FLT8S")

# Compute d8 dirs
dirs <- flowdem::dirs(filled_eps)
terra::writeRaster(dirs, paste0(devpath, "/inst/extdata/dirs.tif"), datatype = "INT1U", overwrite = TRUE)

# Compute d8 flow accumulation
accum <- flowdem::accum(dirs)
terra::writeRaster(accum, paste0(devpath, "/inst/extdata/accum.tif"), overwrite = TRUE)

# Compute watershed of a point
p <- sf::st_read(paste0(devpath, "/inst/extdata/point.gpkg"))
watershed_point <- flowdem::watershed(dirs, p)
terra::writeRaster(watershed_point, paste0(devpath, "/inst/extdata/watershed_point.tif"), datatype = "INT1U", overwrite = TRUE)

# Compute watershed of a line
l <- sf::st_read(paste0(devpath, "/inst/extdata/line.gpkg"))
watershed_line <- flowdem::watershed(dirs, l)
terra::writeRaster(watershed_line, paste0(devpath, "/inst/extdata/watershed_line.tif"), datatype = "INT1U", overwrite = TRUE)

# Compute watershed of a polygon
poly <- sf::st_buffer(p, 500) #500 m buffer
sf::st_write(poly, paste0(devpath, "/inst/extdata/poly.gpkg"))

poly <- sf::st_read(paste0(devpath, "/inst/extdata/poly.gpkg"))
watershed_poly <- flowdem::watershed(dirs, poly)
terra::writeRaster(watershed_poly, paste0(devpath, "/inst/extdata/watershed_poly.tif"), datatype = "INT1U", overwrite = TRUE)

# Compute watershed using a raster
poly_rast <- terra::rasterize(poly, dirs)
terra::writeRaster(poly_rast, paste0(devpath, "/inst/extdata/poly_rast.tif"), datatype = "INT1U", overwrite = TRUE)

poly_rast <- terra::rast(paste0(devpath, "/inst/extdata/poly_rast.tif"))
watershed_poly_rast <- flowdem::watershed(dirs, poly_rast)
terra::writeRaster(watershed_poly_rast, paste0(devpath, "/inst/extdata/watershed_poly_rast.tif"), datatype = "INT1U", overwrite = TRUE)
