% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pf_basins_barnes2014}
\alias{pf_basins_barnes2014}
\title{Improved priority flood with watershed labels (algorithm 5) in:
"Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024"}
\usage{
pf_basins_barnes2014(dem)
}
\arguments{
\item{dem}{The input digital elevation model (DEM)}
}
\value{
List of two rasters: one with the filled input dem and one integer raster with basin labels
}
\description{
Improved priority flood with watershed labels (algorithm 5) in:
"Barnes, R., Lehman, C., Mulla, D., 2014. Priority-flood: An optimal depression-filling and watershed-labeling algorithm for digital elevation models. Computers & Geosciences 62, 117–127. doi:10.1016/j.cageo.2013.04.024"
}
