
# flowdem: Revolutionize your flow routing with digital elevation models in R!

## Overview

Are you tired of struggling with complex algorithms to process and
analyze digital elevation models (DEMs) for accurate water flow
modeling? Look no further! The “flowdem” package is here to simplify
your life. With specialized algorithms, “flowdem” allows you to
delineate watersheds, stream networks, and perform other essential tasks
related to flow routing on DEMs with ease.

## Installation

To enjoy the latest version of “flowdem”, simply use the following
command in R with the “remotes” package:

``` r
remotes::install_github("KennethTM/flowdem")
```

## Key Features

- Efficient DEM processing with algorithms optimized for large grid
  sizes.
- Delineate watersheds (also referred to as drainage basins, catchments,
  upslope areas, or contributing areas) for any location.
- Perform stream and river delineation.
- Accurately route water flows.
- Handle operations in memory, avoiding repeated file writing when
  delineating many watersheds
- And so much more!

## Background

Water flow modeling using DEMs can be a complex process that requires
specialized algorithms to preprocess, analyze, and interpret the data
accurately. “flowdem” simplifies this process by providing you with an
easy-to-use R package that incorporates advanced algorithms from the
[RichDEM library](https://github.com/r-barnes/richdem) by Richard
Barnes. Key algorithms are exposed to R using Rcpp for optimal
performance.

## DEM Processing

After pre-processing your DEM, you can unleash its full potential by
delineating watersheds and streams, routing water flows accurately, and
performing other essential tasks. The pre-processing steps include
dealing with depressions or flat surfaces through *filling* or
*breaching*, enabling correct assessment of flow directions for further
analysis.

## Example Usage

The “flowdem” package works seamlessly with the [“terra”
R-package](https://github.com/rspatial/terra) to manage raster data. If
you encounter any issues, make sure your “terra” package is up-to-date
before reaching out for assistance.

Here”s a quick sample of how “flowdem” can transform your DEM analysis.
The sample data, a DEM (100 m resolution) from northern Zealand,
Denmark, is included in the package.

``` r
library(flowdem);library(terra)

#Read the raster data
dem_path <- system.file("extdata", "arre.tif", package="flowdem")
dem <- rast(dem_path)

plot(dem)
```

![](https://github.com/KennethTM/flowdem/blob/main/man/figures/arre_dem.png)

To accurately determine flow directions, it”s essential to use a DEM
free of depressions and pits. Pre-processing techniques like filling or
breaching can help achieve this:

``` r
#Fill DEM and leave surfaces flat
dem_fill <- fill(dem, epsilon = FALSE)

#Fill DEM and apply a gradient on flat surfaces to ensure flow
#Note: When writing DEMs filled with epsilon to file, 
#set the datatype to FLT8S in terra::writeRaster()
dem_fill_eps <- fill(dem, epsilon = TRUE)

#Fill DEM and delineate coastal drainage basins simultaneously
dem_fill_basins <- fill_basins(dem)

#Breach DEM, that is, resolve depression by "carving" through obstacles
dem_breach <- breach(dem)

#Use fill with epsilon on breached DEM to resolve flats and ensure drainage
dem_breach_fill_eps <- fill(dem_breach, epsilon = TRUE)
```

As demonstrated in this section, the filling and breaching processes
have a significantly different impact on the Digital Elevation Model
(DEM) compared to the raw DEM. To further illustrate this point, we
analyze the differences between the original DEM and the modified
versions.

``` r
dem_breach_diff <- dem - dem_breach
dem_breach_diff[dem_breach_diff == 0] <- NA

dem_fill_diff <- dem - dem_fill_eps
dem_fill_diff[dem_fill_diff == 0] <- NA

par(mfrow = c(1, 2))
plot(dem_breach_diff, main = "Impact of breaching", col=hcl.colors(50, rev = TRUE))
plot(dem_fill_diff, main = "Impact of filling", col=hcl.colors(50, rev = TRUE))
```

![](https://github.com/KennethTM/flowdem/blob/main/man/figures/arre_diff.png)

To improve accuracy in your topographical analysis, properly fill or
breach any sinks and depressions on the DEM. Additionally, process flats
using `epsilon = TRUE` to ensure proper flow direction determination.
Once completed, this information can then be used to effectively outline
the pathways taken by the flowing water.

``` r
#Get flow directions using the filled DEM
dem_dir <- dirs(dem_breach_fill_eps, mode="d8")

#Get flow accumulation
dem_acc <- accum(dem_dir, mode="d8")

#Flow accumulation can be used for stream delineation, 
#e.g using a threshold of 100 contributing cells (100*100*100 = 1 km2)
dem_streams <- dem_acc > 100

plot(dem_streams, col=c("grey", "dodgerblue"), legend=FALSE)
```

![](https://github.com/KennethTM/flowdem/blob/main/man/figures/arre_streams.png)

The flow directions are based on the “d8” routing scheme, which is
currently the only mode available. Utilizing the flow direction grid can
help with watershed delineation for various applications. You can
achieve this by providing a target area in either raster or vector
format from the “terra” package or vector format from the “sf” package.
In this example, we will demonstrate how to delineate the watershed for
Denmark’s largest lake (Lake Arre).:

``` r
#Read the vector data using terra of sf packages
lake_path <- system.file("extdata", "arre.sqlite", package="flowdem")

#Using sf
# lake <- st_read(lake_path) 

#Using terra
lake <- vect(lake_path)

#Using terra raster
# lake <- rasterize(lake, dem_dir)

lake_watershed <- watershed(dem_dir, lake)

plot(dem)
plot(lake_watershed, col="dodgerblue", add=TRUE, legend=FALSE, alpha=0.5)
plot(lake, col="coral", add=TRUE, border="coral")
```

![](https://github.com/KennethTM/flowdem/blob/main/man/figures/arre_watershed.png)

The library provides support for creating watershed delineations
targeting points, lines or polygons in terra” or “sf” packages. It also
supports the use of rasters in “terra”. Additionally, this library
allows for nested watershed delineation with each watershed assigned a
unique ID.

## Contributors

[Cyril Mory](https://github.com/cyrilmory)
