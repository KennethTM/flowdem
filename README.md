
# flowdem package

## R-package for flow routing on digital elevation models

### Installation

Use the ‘remotes’ package to install the package from Github:

``` r
remotes::install_github('KennethTM/flowdem')
```

### Background

Water flow can be modelled using a digital elevation model (DEM) often
found in the form of a raster or grid where each grid cell has an
elevation value. Using DEM for flow routing requires specialized
algorithms for preprocessing, routing and analysis that should be
suitable for large grids considering the increasing resolution of DEMs
collected today.

This R-package provides the tools for processing DEMs. Code for the main
algorithms are ported from the RichDEM library by Richard Barnes which
has been implemented in this package using Rcpp. For more information
and algorithms see the [RichDEM Python and C++
libraries](https://github.com/r-barnes/richdem).

### DEM processing

After processing the DEM, several useful things can be done including:
delineation of a watershed (also termed drainage basin, catchment,
upslope or contributing area) for any location, stream and river
delineation, routing water flow, and much more. The pre-processing often
involve steps that deal with depression and/or flat surfaces by either
*filling* or *breaching* them to enable correct assessment of *flow
directions*. Flow directions are then used for watershed and stream
delineation etc.

### Example useage

The library uses the [‘terra’
R-package](https://github.com/rspatial/terra) for managing raster data.
If you encounter any issues, first try to update ‘terra’.

Here, a sample DEM from Denmark (100 m resolution) included in the
package is used to showcase functions of the library:

``` r
library(flowdem);library(terra)

#Read the raster data
dem_path <- system.file("extdata", "arre.tif", package="flowdem")
dem <- rast(dem_path)

plot(dem)
```

![](https://github.com/KennethTM/flowdem/blob/main/man/figures/arre_dem.png)

Determining flow directions requires a DEM without depression and pits.
Filling and/or breaching can be used to pre-process the DEM:

``` r
#Fill DEM and leave surfaces flat
dem_fill <- fill(dem, epsilon = FALSE)

#Fill DEM and apply a gradient on flat surfaces to ensure flow
#Note: When writing DEMs filled with epsilon to file, set the datatype to FLT8S in terra::writeRaster()
dem_fill_eps <- fill(dem, epsilon = TRUE)

#Fill DEM and delineate coastal drainage basins simultaneously
dem_fill_basins <- fill_basins(dem)

#Breach DEM, that is, resolve depression by 'carving' through obstacles
dem_breach <- breach(dem)

#Use fill with epsilon on breached DEM to resolve flats and ensure drainage
dem_breach_fill_eps <- fill(dem_breach, epsilon = TRUE)
```

The filling and breaching operations affects the DEM very differently as
shown here where differences between the raw DEM and the filled and
breached DEMs in a small sub-region:

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

After filling and/or breaching and resolving flats on the DEM, flow
directions can be calculated and used to delineate flow paths:

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

The flow directions are based on the ‘d8’ routing scheme which is the
only mode supported now. Using the flow direction grid to delineate
watersheds are useful for many applications. This can be done by
supplying a target area in the form of raster or data in a spatial
vector format through the ‘sp’ or ‘sf’ packages. In this example the
watershed for the largest lake in Denmark (Lake Arre) is delineated:

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

The library support watershed delineation to a target point, line or
polygon used in the ‘terra’ or ‘sf’ packages. A raster can also be
supplied. Furthermore, the library supports delineation of nested
watersheds where each watershed is assigned a unique id.

### Contributors

[Cyril Mory](https://github.com/cyrilmory)
