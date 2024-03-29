#Functions for dealing with depressions in digital elevation models

#' Remove depressions by filling
#'
#' Remove depressions from a digital elevation model by filling it inwards from the edges using the Priority-Flood algorithm. See details to be aware of when write the results to file.
#'
#' WARNING: The correction added to ensure flow on flat surfaces is REALLY small (< 1e-8)
#' which is not a problem as long as filled_eps remains in memory (terra always uses FLT8S for cell values
#' while in memory), but may be lost during writing to disk (terra uses datatype set in terraOptions()
#' for cell values when writing to disk, default is FLT4S). When saving DEMs filled with epsilon to disk,
#' set the datatype to FLT8S either in terraOptions() or in terra::writeRaster().
#'
#' @md
#' @param dem terra::SpatRaster object containing the digital elevation model.
#' @param epsilon TRUE (default) or FALSE. If TRUE, cell elevations in depressions are be increased to ensure drainage. If FALSE, filled depressions are left as flat surfaces.
#' @return dem_fill terra::SpatRaster object containing the digital elevation model with depressions removed.
#' @export fill 
#' @export
fill <- function(dem, epsilon = TRUE){
  
  if(!inherits(dem, "SpatRaster")){
    stop("Input must be a SpatRaster object from the terra package")
  }

  dem_mat <- terra::as.matrix(dem, wide=TRUE)
  class(dem_mat) <- "numeric" #explicit conversion to double to avoid issues when integer matrix are passed
  dem_mat[is.na(dem_mat)] <- -9999

  if(epsilon){
    pf_eps_barnes2014(dem_mat)
  }else{
    pf_barnes2014(dem_mat)
  }

  dem_mat[dem_mat == -9999] <- NA
  terra::values(dem) <- dem_mat

  return(dem)
}

#' Remove depressions by filling and delineate drainage basins
#' 
#' Remove depressions from a digital elevation model by filling it inwards from the edges using the Priority-Flood algorithm, and delineate drainage basins simultaneously.
#' 
#' @md
#' @param dem terra::SpatRaster object containing the digital elevation model.
#' @return dem_fill_basins terra::SpatRaster object with two layers: one with the filled input dem and one integer raster with basin labels
#' @export fill_basins 
#' @export
fill_basins <- function(dem){

  if(!inherits(dem, "SpatRaster")){
    stop("Input must be a SpatRaster object from the terra package")
  }

  dem_mat <- terra::as.matrix(dem, wide=TRUE)
  class(dem_mat) <- "numeric"
  dem_mat[is.na(dem_mat)] <- -9999

  mat_list <- pf_basins_barnes2014(dem_mat)

  mat_list$dem[mat_list$dem == -9999] <- NA
  terra::values(dem) <- mat_list$dem

  mat_list$labels[mat_list$labels == 0] <- NA
  dem_basins <- dem
  terra::values(dem_basins) <- mat_list$labels

  dem_fill_basins <- terra::rast(list(dem, dem_basins))
  names(dem_fill_basins) <- c("dem", "basins")

  return(dem_fill_basins)

}

#' Remove depressions by breaching
#'
#' Remove depressions from digital elevation models by breaching depressions
#' 
#' @md
#' @param dem terra::SpatRaster object containing the digital elevation model.
#' @return dem_breach terra::SpatRaster object with depressions filled.
#' @export breach 
#' @export
breach <- function(dem){
  
  if(!inherits(dem, "SpatRaster")){
    stop("Input must be a SpatRaster object from the terra package")
  }
  
  dem_mat <- terra::as.matrix(dem, wide=TRUE)
  class(dem_mat) <- "numeric"
  dem_mat[is.na(dem_mat)] <- -9999
  
  comp_breach_lindsay2016(dem_mat)
  
  dem_mat[dem_mat == -9999] <- NA
  terra::values(dem) <- dem_mat
  
  return(dem)
  
}
