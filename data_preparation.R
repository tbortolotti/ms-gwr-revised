# Working directory (change)
setwd("/home/giovanni/Scrivania/teresa/ms-gwr-revised")

rm(list=ls())
graphics.off()
cat("\014")

## DATA PREPARATION
## This script is used to projects event and site coordinates on a UTM system, and prepares the dataset
## that is used for the calibration of the Geographically weighted model.
## It loads data from files:
## - "data_dir/italian_data_pga.RData"
## - "data_dir/confini_ut33.shp"
## and saves the pre-processed data in
## two files:
## - "data_dir/utm_coordinates.RData"
## - "data_dir/regressors.RData"

# Load data
dataset = readRDS("data_dir/italian_data_pga.RData")

# write.csv(dataset, "data_dir/italian_data_pga.csv")

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

# Remove Sortino
SRT.lat = 37.16657
SRT.long = 15.05421

SRT.idx = 4632
dataset = dataset[-SRT.idx,]

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

# Midpoint
utm_md = (utm_ev + utm_st)/2

# Build spatial points data frame for UTM33 coordinates
utm_ev_sp = SpatialPointsDataFrame(utm_ev, dataset[,1:6])
utm_st_sp = SpatialPointsDataFrame(utm_st, dataset[,1:6])
utm_md_sp = SpatialPointsDataFrame(utm_md, dataset[,1:6])

# Transform shapefile into spatial dataframe
shape_utm_spatial = as_Spatial(shape_utm_no_sardinia)
grid_utm = makegrid(shape_utm_spatial, cellsize = 10000) # cellsize in map units!
grid_utm = SpatialPoints(grid_utm, proj4string = CRS(proj4string(shape_utm_spatial)))
grid_inside_utm = grid_utm[shape_utm_spatial,]
coords_utm = grid_inside_utm@coords
coords_df_utm = as.data.frame(coords_utm)

save(utm_md_sp, utm_ev_sp, utm_st_sp, coords_utm, file="data_dir/utm_coordinates.RData")
save(utm_md, utm_ev, utm_st, shape_utm_no_lamp, file="data_dir/spacial_info_plots.RData")


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

save(b1,b2,c1,c2,c3,f1,f2,k,mh,mref,h,y, file="data_dir/regressors.RData")

