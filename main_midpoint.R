# Working directory (change)
setwd("C:/Users/Teresa Bortolotti/Documents/R/ms-gwr-revised")

source('package_install.R')

library(sf)
library(raster)
library(ggplot2)
library(rgdal) 
library(GWmodel)
library(cowplot)
library(geosphere)
library(psych)
library(pracma)
library(reshape2)
library(plot3D)
library(progress)
library(roahd)
library(ggplot2)

rm(list=ls())
graphics.off()
cat("\014")

## Load -------------------------------------------------
# Load data
dataset = readRDS("data_dir/italian_data_pga.RData")

# Load functions
#source("functions/functions.R")
source("functions/gcv_mei_only_one.R")
source("functions/SEC_only_calibration.R")
source("functions/ESC_only_calibration.R")
source("functions/SEC_only_constant_intercept_calibration.R")
source("functions/SEC_no_intercept_calibration.R")
source("functions/SEC_grid_creation.R")
source("functions/gauss_kernel.R")
source("functions/stationarity_check.R")
source("functions/significance_check.R")

# shapefiles
shape_utm = st_read('data_dir/confini_ut33.shp')

## Shapefiles -------------------------------------------

# Removal of unused regions
shape_utm_no_lamp = shape_utm
shape_utm_no_lamp$geometry[[1]][[12]] = NULL

shape_utm_no_sardinia = shape_utm
shape_utm_no_sardinia$geometry[[1]][[3]] = NULL
for (i in 4:10){
  shape_utm_no_sardinia$geometry[[1]][[4]] = NULL
}
for (i in 12:57){
  shape_utm_no_sardinia$geometry[[1]][[5]] = NULL
}

# Projection of event and site coordinates to work with UTM33
latitude_ev = dataset$ev_latitude
longitude_ev = dataset$ev_longitude
long_lat_ev = cbind(longitude_ev, latitude_ev)
utm_ev = project(long_lat_ev, "+proj=utm +zone=33 ellps=WGS84") 
utm_ev = as.data.frame(utm_ev)
long_lat_ev = as.data.frame(long_lat_ev)

latitude_st = dataset$st_latitude
longitude_st = dataset$st_longitude
long_lat_st = cbind(longitude_st, latitude_st)
utm_st = project(long_lat_st, "+proj=utm +zone=33 ellps=WGS84") 
utm_st = as.data.frame(utm_st)
long_lat_st = as.data.frame(long_lat_st)

# Build spatial points data frame for UTM33 coordinates
utm_ev_sp = SpatialPointsDataFrame(utm_ev, dataset[,1:6])
utm_st_sp = SpatialPointsDataFrame(utm_st, dataset[,1:6])

# Transform shapefile into spatial dataframe
shape_utm_spatial = as_Spatial(shape_utm_no_sardinia)
grid_utm = makegrid(shape_utm_spatial, cellsize = 10000) # cellsize in map units!
grid_utm = SpatialPoints(grid_utm, proj4string = CRS(proj4string(shape_utm_spatial)))
grid_inside_utm = grid_utm[shape_utm_spatial,]
coords_utm = grid_inside_utm@coords
coords_df_utm = as.data.frame(coords_utm)

## Build the data -------------------------------------
mh = 5.5
mref = 5.324
h = 6.924
attach(dataset)
b1 = (mag-mh)*(mag<=mh)
b2 = (mag-mh)*(mag>mh)
c1 = (mag-mref)*log10(sqrt(JB_complete^2+h^2))
c2 = log10(sqrt(JB_complete^2+h^2))
c3 = sqrt(JB_complete^2+h^2)
f1 = as.numeric(fm_type_code == "SS")
f2 = as.numeric(fm_type_code == "TF")
k = log10(vs30/800)*(vs30<=1500)+log10(1500/800)*(vs30>1500)
y = log10(rotD50_pga)
detach(dataset)

Xc = cbind(b1,b2,f1,f2,c1)
Xe = cbind(c2,c3)
Xs = k
