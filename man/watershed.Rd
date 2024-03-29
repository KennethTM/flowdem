% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/d8_functions.R
\name{watershed}
\alias{watershed}
\title{Delineate watersheds}
\usage{
watershed(dirs, target, nested = FALSE, mode = "d8")
}
\arguments{
\item{dirs}{terra::SpatRaster object containing d8 flow direction.}

\item{target}{terra::SpatRaster, terra::SpatVector or sf::sf object.}

\item{nested}{TRUE or FALSE (default). Indicates whether the output watersheds should be nested (a watershed for each pour-point).}

\item{mode}{Only 'd8' supported for now.}
}
\value{
watershed terra::SpatRaster object delineated watershed.
}
\description{
Delineate watersheds (equivalent terms: catchment area, upslope area, contributing area, basin) from a flow direction raster and target area given as vector or raster object. For rasters, all non-NA cells are used as pour points.
}
