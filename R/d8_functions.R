#Functions using the d8 flow routing model

#' Determine flow directions on digital elevation models
#' 
#' @md
#' @param dem terra::SpatRaster object containing the digital elevation model.
#' @param mode Only 'd8' supported for now.
#' @return dirs terra::SpatRaster object with flow directions.
#' @export dirs 
#' @export
dirs <- function(dem, mode = "d8"){
  
  if(!inherits(dem, "SpatRaster")){
    stop("Input must be a SpatRaster object from the terra package")
  }
  
  if(mode != "d8"){
    stop("Only the 'deterministic eight' (d8) flow model is supported for now.")
  }
  
  dem_mat <- terra::as.matrix(dem, wide=TRUE)
  class(dem_mat) <- "numeric"
  dem_mat[is.na(dem_mat)] <- -9999
  
  dirs_mat <- d8_flow_directions(dem_mat)
  
  dirs_mat[dirs_mat == 0] <- NA

  dirs <- dem
  terra::values(dirs) <- dirs_mat
  
  return(dirs)
  
}

#' Determine flow accumulation on digital elevation models
#' 
#' @md
#' @param dirs terra::SpatRaster object with flow directions.
#' @param mode Only 'd8' supported for now.
#' @return accum terra::SpatRaster object with flow accumulation.
#' @export accum 
#' @export
accum <- function(dirs, mode = "d8"){
  
  if(!inherits(dirs, "SpatRaster")){
    stop("Input must be a SpatRaster object from the terra package")
  }
  
  if(mode != "d8"){
    stop("Only the 'deterministic eight' (d8) flow model is supported for now.")
  }
  
  mm <- terra::minmax(dirs, compute = TRUE)
  input_min <- mm[1]
  input_max <- mm[2]

  if(!is.integer(dirs[]) | input_min < 1 | input_max > 8){
    stop("Input must have d8 flow directions 
          encoded as integers between 1 and 8
          using the following flow coordinate 
          system (0 is the focal cell):
          234
          105
          876")
  }
  
  dirs_mat <- terra::as.matrix(dirs, wide=TRUE)
  dirs_mat[is.na(dirs_mat)] <- 0
  
  acc_mat <- d8_flow_accum(dirs_mat)
  
  acc_mat[acc_mat == -1] <- NA
  accum <- dirs
  terra::values(accum) <- acc_mat
  
  return(accum)
  
}


#' Delineate watersheds (equivalent terms: catchment area, upslope area, contributing area, basin) from a flow direction raster and target area given as vector or raster object.
#' 
#' @md
#' @param dirs RasterLayer object containing d8 flow direction.
#' @param target RasterLayer, Spatial* or sf object.
#' @param nested TRUE or FALSE (default). Indicates whether the output watersheds should be nested (a watershed for each pour-point).
#' @param mode Only 'd8' supported for now.
#' @return watershed RasterLayer object delineated watershed.
#' @export watershed 
#' @export
watershed <- function(dirs, target, nested = FALSE, mode = "d8"){
  
  if(!any(inherits(target, "RasterLayer"), !inherits(target, "sf"), !inherits(target, "Spatial"))){
    stop("Input must be a RasterLayer, Spatial* or sf object")
  }
  
  if(mode != "d8"){
    stop("Only the 'deterministic eight' (d8) flow model is supported for now.")
  }
  
  if(dirs@data@haveminmax){
    input_min <- dirs@data@min
    input_max <- dirs@data@max
  }else{
    input_min <- min(dirs[], na.rm = TRUE)
    input_max <- max(dirs[], na.rm = TRUE)
  }
  
  if(!is.integer(dirs[]) | input_min < 1 | input_max > 8){
    stop("Input must have d8 flow directions 
          encoded as integers between 1 and 8
          using the following flow coordinate 
          system (0 is the focal cell):
          234
          105
          876")
  }
  
  if(inherits(target, "sf") | inherits(target, "Spatial")){
    
    # Special case when target is a single linestring
    if((nrow(target) == 1) & ((sf::st_geometry_type(target)[1]) == "LINESTRING")){
      target_cells <- raster::extract(dirs, target, cellnumbers = TRUE)
      target_xy <- raster::rowColFromCell(dirs, target_cells[[1]][,1])
      target_xy <- cbind(target_xy, rep(1, nrow(target_xy)))
    } else {
      # General case (crashes when target is a single line)
      cell_df <- raster::extract(dirs, target, cellnumbers=TRUE, df = TRUE)
      target_xy <- raster::rowColFromCell(dirs, cell_df$cell)
      target_xy <- cbind(target_xy, cell_df$ID)
    }
	
  }else if(inherits(target, "RasterLayer")){
    
    if(!raster::compareRaster(dirs, target)){
      stop("Input dirs and target rasters must match")
    }
    
    target_cells <- which(!is.na(target[]))
    target_xy <- raster::rowColFromCell(dirs, target_cells)
    target_xy <- cbind(target_xy, target[target_cells])
    
  }else{
    stop("Something went wrong")
  }
  
  if(nrow(target_xy) == 0){
    stop("No outlets found. 
          Do the function arguments dirs and target intersect geographically? Do the projections match each other?")
  }
  
  dirs_mat <- raster::as.matrix(dirs)
  dirs_mat[is.na(dirs_mat)] <- 0
  
  watershed_mat <- d8_watershed_nested(dirs_mat, target_xy, nested = nested)
  watershed_mat[watershed_mat == 0] <- NA
  watershed <- raster::raster(watershed_mat, template = dirs)
  raster::dataType(watershed) <- "INT4U"
  
  return(watershed)
  
}
