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
  
  return(raster_copy)
  
}

#' Remove depressions from digital elevation models by filling using the Priority-Flood algorithm and get watershed labels
#' 
#' @md
#' @param raster RasterLayer object.
#' @return stack RasterStack object two layers: one with the filled input dem and one integer raster with watershed labels
#' @export dem_fill_watersheds 
#' @export
dem_fill_watersheds <- function(raster){
  
  if(!inherits(raster, "RasterLayer")){
    stop("Input must be a RasterLayer object from the raster package")
  }
  
  #copy raster to avoid inplace modification or use Rcpp clone()
  raster_copy <- raster
  
  raster_copy[is.na(raster_copy)] <- -9999
  
  dem <- raster::as.matrix(raster_copy)
  
  dem_list <- pf_watersheds_barnes2014(dem)
  
  raster_copy[] <- dem_list$dem
  raster_copy[raster_copy == -9999] <- NA
  
  raster_copy_2 <- raster
  raster_copy_2[] <- dem_list$labels
  raster_copy_2[raster_copy_2 == 0] <- NA
  
  raster_stack <- stack(raster_copy, raster_copy_2)
  names(raster_stack) <- c("dem", "watersheds")
  
  return(raster_stack)
  
}
