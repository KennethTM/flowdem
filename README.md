# flowdem package

## R-package for flow routing on digital elevation models
## *Under ongoing development and testing*

### Installation

Use the 'remotes' package to install the package from Github:

```r
remotes::install_github('KennethTM/flowdem')
```

### Background

Water flow can be modelled using a digital elevation model (DEM) often found in the form of a raster or grid where each grid cell has an elevation value. Using DEM for flow routing requires specialized algorithms for preprocessing, routing and analysis that should be suitable for large grids considering the increasing resolution of DEMs collected today. 

This R-package provides the tools for processing DEMs. Code for the main algorithms are ported from the RichDEM library by Richard Barnes which has been implemented in this package using Rcpp. For more information and algorithms see the [RichDEM Python and C++ libraries](https://github.com/r-barnes/richdem). 

### DEM processing

After processing the DEM, several useful things can be done including: delineation of a watershed (also termed drainage basin, catchment, upslope or contributing area) for any location, stream and river delineation, routing water flow, and much more. The pre-processing often involve steps that deal with depression and/or flat surfaces by either *filling* or *breaching* them to enable correct assessment of *flow directions*. Flow directions are then used for watershed and stream delineation etc. 

### Example useage

The library uses the 'terra' R-package for managing raster data. Here, a sample DEM from Denmark (10 m resolution) is used to showcase functions of the library:

```r
library(flowdem);library(terra)

#Read in the raster data
dk_dem <- rast("test/dk_dem_sub.tif")

#Fill DEM and leave surfaces flat
dk_dem_fill <- fill(dk_dem, epsilon = FALSE)

#Fill DEM and apply a gradient on flat surfaces to ensure flow
dk_dem_fill_eps <- fill(dk_dem, epsilon = TRUE)

#Fill DEM and delineate coastal drainage basins simultenously
dk_dem_fill_basins <- fill_basins(dk_dem)

#Breach DEM, that is, resolve depression by 'carving' through obstacles
dk_dem_breach <- breach(dk_dem)

#Use fill with epsilon on breached DEM to resolve flats and ensure drainage
dk_dem_breach_fill_eps <- fill(dk_dem_breach, epsilon = TRUE)
```

The filling and breaching operations affects the DEM very differently as shown here where differences between the raw DEM and the filled and breached DEMs in a small sub-region:

![](https://github.com/KennethTM/flowdem/blob/main/test/fill_vs_breach.png)

After filling and/or breaching and resolving flats on the DEM, flow directions can be calculated and used to delineate flow paths for example:

```r
#Get flow directions using the filled DEM
dk_dem_dir <- dirs(dk_dem_breach_fill_eps, mode="d8")

#Get flow accumulation which can be used for stream delineation, e.g using a threshold (dk_dem_acc > 1000)
dk_dem_acc <- accum(dk_dem_breach_fill_eps, mode="d8")
```

The flow directions are based on the 'd8' routing scheme which is the only mode supported now. Using the flow direction grid to delineate watersheds are useful for many applications. This can be done by supplying a target area in the form of raster or data in a spatial vector format through the 'sp' or 'sf' packages. In this example the watershed for the largest lake in Denmark is delineated:

```r
library(sf)

lake <- st_read("test/lake_arre.sqlite")

lake_watershed <- watershed(dk_dem_breach_fill_eps, lake)
```

![](https://github.com/KennethTM/flowdem/blob/main/test/watershed_delin.png)

The library support watershed delineation to a target point, line or polygon used in the 'sp' or 'sf' packages. A raster can also be supplied. Furthermore, the library supports delineation of nested watersheds where each watershed is assigned a unique id.

### Future development

Goals for future development:

* Add function for flat resolution
* Add latest algorithms of the 'priority flood family'
* Improve error messages
* Add weight option for flow accumulation
* Add support for dInf flow routing
