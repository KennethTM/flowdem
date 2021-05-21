#' Remove depressions from digital elevation models by filling using the Priority-Flood algorithm
#' 
#' @md
#' @param raster RasterLayer object.
#' @param epsilon TRUE (default) or FALSE. Indicates whether cell elevations in depressions should be increased to ensure drainage.
#' @return raster_fill RasterLayer object with depressions filled.
#' @export dem_fill 
#' @export
dem_fill <- function(raster, epsilon = TRUE){
  
  if(!inherits(raster, "RasterLayer")){
    stop("Input must be a RasterLayer object from the raster package")
  }
  
  #copy raster to avoid inplace modification or use Rcpp clone()
  raster_copy <- raster
  
  raster_copy[is.na(raster_copy)] <- -9999
  
  dem <- raster::as.matrix(raster_copy)
  
  if(epsilon){
    dem_fill <- pf_eps_barnes2014(dem)
  }else{
    dem_fill <- pf_barnes2014(dem)
  }
  
  raster_copy[] <- dem_fill
  raster_copy[raster_copy == -9999] <- NA
  raster::dataType(raster_copy) <- "FLT8S"
  
  return(raster_copy)
  
}

#' Remove depressions from digital elevation models by filling using the Priority-Flood algorithm and get basin labels
#' 
#' @md
#' @param raster RasterLayer object.
#' @return stack RasterStack object two layers: one with the filled input dem and one integer raster with basin labels
#' @export dem_fill_basins 
#' @export
dem_fill_basins <- function(raster){
  
  if(!inherits(raster, "RasterLayer")){
    stop("Input must be a RasterLayer object from the raster package")
  }
  
  #copy raster to avoid inplace modification or use Rcpp clone()
  raster_copy <- raster
  
  raster_copy[is.na(raster_copy)] <- -9999
  
  dem <- raster::as.matrix(raster_copy)
  
  dem_list <- pf_basins_barnes2014(dem)
  
  raster_copy[] <- dem_list$dem
  raster_copy[raster_copy == -9999] <- NA
  raster::dataType(raster_copy) <- "FLT8S"
  
  raster_copy_2 <- raster
  raster_copy_2[] <- dem_list$labels
  raster_copy_2[raster_copy_2 == 0] <- NA
  raster::dataType(raster_copy_2) <- "INT2U"
  
  raster_stack <- stack(raster_copy, raster_copy_2)
  names(raster_stack) <- c("dem", "basins")
  
  return(raster_stack)
  
}

#' Remove depressions from digital elevation models by breaching depressions
#' 
#' @md
#' @param raster RasterLayer object.
#' @return raster_fill RasterLayer object with depressions filled.
#' @export dem_breach 
#' @export
dem_breach <- function(raster){
  
  if(!inherits(raster, "RasterLayer")){
    stop("Input must be a RasterLayer object from the raster package")
  }
  
  #copy raster to avoid inplace modification or use Rcpp clone()
  raster_copy <- raster
  
  raster_copy[is.na(raster_copy)] <- -9999
  
  dem <- raster::as.matrix(raster_copy)
  
  dem_breach <- comp_breach_lindsay2016(dem)
  
  raster_copy[] <- dem_breach
  raster_copy[raster_copy == -9999] <- NA
  raster::dataType(raster_copy) <- "FLT8S"
  
  return(raster_copy)
  
}

#' Determine flow directions on digital elevation models
#' 
#' @md
#' @param raster RasterLayer object.
#' @param mode Only 'd8' supported for now.
#' @return raster_direction RasterLayer object with flow directions.
#' @export dem_direction 
#' @export
dem_direction <- function(raster, mode = "d8"){
  
  if(!inherits(raster, "RasterLayer")){
    stop("Input must be a RasterLayer object from the raster package")
  }
  
  if(mode != "d8"){
    stop("Only the 'deterministic eight' (d8) flow model is supported for now.")
  }
  
  #copy raster to avoid inplace modification or use Rcpp clone()
  raster_copy <- raster
  
  raster_copy[is.na(raster_copy)] <- -9999
  
  dem <- raster::as.matrix(raster_copy)
  
  dem_dirs <- d8_flow_directions(dem)
  
  raster_copy[] <- dem_dirs
  raster_copy[raster_copy == 0] <- NA
  raster::dataType(raster_copy) <- "INT1U"
  
  return(raster_copy)
  
}

#' Determine flow accumulation on digital elevation models
#' 
#' @md
#' @param raster RasterLayer object.
#' @param mode Only 'd8' supported for now.
#' @return raster_accumulate RasterLayer object with flow accumulation.
#' @export dem_accumulation 
#' @export
dem_accumulation <- function(raster, mode = "d8"){
  
  if(!inherits(raster, "RasterLayer")){
    stop("Input must be a RasterLayer object from the raster package")
  }
  
  if(mode != "d8"){
    stop("Only the 'deterministic eight' (d8) flow model is supported for now.")
  }
  
  #copy raster to avoid inplace modification or use Rcpp clone()
  raster_copy <- raster
  
  raster_copy[is.na(raster_copy)] <- -9999
  
  dem <- raster::as.matrix(raster_copy)
  
  dem_acc <- d8_flow_accum(dem)
  
  raster_copy[] <- dem_acc
  raster_copy[raster_copy == -1] <- NA
  raster_copy[is.na(raster[])] <- NA
  raster::dataType(raster_copy) <- "FLT8S"
  
  return(raster_copy)
  
}


#' Delineate watersheds (same as catchment, upslope area, contributing area) from a flow direction raster and target area given as point, polygon or raster object.
#' 
#' @md
#' @param direction RasterLayer object containing d8 flow direction.
#' @param target RasterLayer, or sf POINT, POLYGON or MULTIPOLYGON object.
#' @param vector_result TRUE or FALSE (default). Indicates whether the output watershed should be converted to a vector geometry.
#' @param nested TRUE or FALSE (default). Indicates whether the output watersheds should be nested (a watershed for each pour-point).
#' @return raster_watershed RasterLayer object delineated watershed.
#' @export dem_watershed 
#' @export
dem_watershed <- function(direction, target, vector_out = FALSE, nested = FALSE){
  
  if(inherits(target, "sf")){
    target_sp <- as(target, "Spatial")
    if(st_is(target, "POINT")){
      target_cells <- cellFromXY(direction, target_sp)
      target_xy <- rowColFromCell(direction, unlist(target_cells))
    }else if(st_is(target, "POLYGON") | st_is(target, "MULTIPOLYGON")){
      target_cells <- cellFromPolygon(direction, target_sp)
      target_xy <- rowColFromCell(direction, unlist(target_cells))
    }else{
      stop("sf type error")
    }
  }else if (inherits(target, "RasterLayer")){
    target_cells <- which(target[] > 0)
    target_xy <- rowColFromCell(direction, unlist(target_cells))
  }else{
    stop("target type error")
  }
  
  direction[is.na(direction)] <- 0
  direction_mat <- raster::as.matrix(direction)
  
  watershed_mat <- d8_watershed(direction_mat, target_xy)
  
  if(nested){
    target_xy <- cbind(target_xy, 1:nrow(target_xy))
    watershed_mat <- d8_watershed_nested(direction_mat, target_xy)
  }
  
  watershed_result <- raster::raster(watershed_mat, template = direction)
  
  if(vector_out){
    watershed_result <- raster::rasterToPolygons(watershed_result, fun=function(x){x>0}, dissolve=TRUE, n=4) #8-connect?
    watershed_result <- sf::st_as_sf(watershed_result)
  }
  
  return(watershed_result)
  
}
