% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{d8_watershed_nested}
\alias{d8_watershed_nested}
\title{Function for d8 watersheds to a target area identified by row-col indexes
Potentially with labeling of nested watersheds}
\usage{
d8_watershed_nested(flowdirs, target_rc, nested)
}
\arguments{
\item{flowdirs}{The d8 pointer flow direction raster}

\item{target_rc}{The outlet}

\item{nested}{Boolean}
}
\value{
a flow accumulation raster
}
\description{
Function for d8 watersheds to a target area identified by row-col indexes
Potentially with labeling of nested watersheds
}
