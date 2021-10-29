
# Working directory (change)
setwd("C:/Users/Teresa Bortolotti/Documents/R/ms-gwr-reviewed")

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

rm(list=ls())
graphics.off()
cat("\014")

## Load -------------------------------------------------
dataset = readRDS("italian_data_pga.RData")
source("functions.R")

# shapefiles
shape_utm = st_read('confini_ut33.shp')

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

## Selection of the optimal bandwidth -------------------------------
bwe_tot = c(10000, 25000, 50000, 75000, 100000)
bws_tot = c(10000, 25000, 50000, 75000, 100000)
n_bwe = length(bwe_tot)
n_bws = length(bws_tot)
bw_table = matrix(0, n_bwe, n_bws)
row.names(bw_table) = bwe_tot 
colnames(bw_table) = bws_tot

for (i in 1:n_bwe){ #i=1
  for (j in 1:n_bws){ # j=1
    bandwidth = gcv_mei_only_one(bwe       = i,
                                 bws       = j,
                                 func      = SEC_only_calibration,
                                 Xc        = Xc,
                                 Xe        = Xe,
                                 Xs        = Xs,
                                 y         = y,
                                 intercept = "c",
                                 utm_ev_sp = coordinates(utm_ev_sp),
                                 utm_st_sp = coordinates(utm_st_sp))
    bw_table[i,j] = bandwidth$gcv
  }
}



bw_best = which(bw_table == max(bw_table), arr.ind = T)
bwe = bwe_tot[bw_best[1]]
bws = bws_tot[bw_best[2]]

bwe = 25000
bws = 75000

## PERMUTATION TESTS ---------------------------
## Comparison of a constant coefficients model with a spatially varying coefficients model -------------
# Joint test for the stationarity of the coefficients
# H0: all coefficients are constant
# H1: at least one coefficient is spatially varying
ols = lm(y ~ Xc + Xe + Xs)

(Start.Time <- Sys.time())
only_intercept = SEC_only_constant_intercept_calibration(Xe        = Xe,
                                                         Xs        = cbind(Xc, Xs),
                                                         y         = y,
                                                         bwe       = bwe,
                                                         bws       = bws,
                                                         utm_ev_sp = coordinates(utm_ev_sp),
                                                         utm_st_sp = coordinates(utm_st_sp))

End.Time <- Sys.time()
## Run-time:
round(End.Time - Start.Time, 2)

save(only_intercept, file="res/only_intercept.RData")
load("res/only_intercept.RData")

n_sample = length(y)
#compute R(H0)
X = cbind(rep(1,n_sample), Xc, Xe, Xs)
I = diag(1,n_sample)
Hols = X%*%(solve(t(X)%*%X))%*%t(X)
RH0 = t(I-Hols)%*%(I-Hols)
epsilon = (I-Hols)%*%y
#compute R(H1)
B = only_intercept$B
#load("lucares/B.RData")
Xcc = rep(1,n_sample)
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
save(RH1, file="RH1_only_intercept_rotD50pga.RData")
#use this RH1 for all the following stationarity tests

#compute T
T0 = (t(y) %*% (RH0-RH1) %*% y) / (t(y) %*% RH1 %*% y)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
print("Permutations")
pb = progress_bar$new(total=n_perm, format = "  computing [:bar] :percent eta: :eta")
pb$tick(0)
for (i in 1:n_perm){
  eps_star = sample(epsilon)
  y_star = Hols %*% y + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
  pb$tick()
}
p_ols = sum(t_stat>as.numeric(T0))/n_perm
p_ols

## Check for the stationarity of coefficients, one at a time ---------------------------
source("functions/stationarity_check.R")
coef_to_check = "b1" #change it among: {"b1","b2","c1","c2","c3","f1","f2","k"}
regs = list(b1 = b1,
            b2 = b2,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            f1 = f1,
            f2 = f2,
            k  = k)
p_one_at_a_time = stationarity_check(coef_to_check = coef_to_check,
                                     regs          = regs,
                                     y             = y,
                                     bwe           = bwe,
                                     bws           = bws,
                                     utm_ev_sp     = utm_ev_sp,
                                     utm_st_sp     = utm_st_sp)
p_one_at_a_time

## Check jointly the stationarity of regressors that are found to be stationary in the
## previous tests -----------------------------------------------------------------------
# H0: b1, b2, f1, f2, c1, other than the intercept, are constant
# H1: at least one is spatially varying
Xc = cbind(b1,b2,f1,f2,c1)
Xe = cbind(c2,c3)
Xs = k
ols = lm(y ~ Xc + Xe + Xs)
sec = SEC_only_calibration(Xc, Xe, Xs, y, "c", bwe, bws, coordinates(utm_ev_sp), coordinates(utm_st_sp))
#compute R(H0)
X = cbind(rep(1,n_sample), Xc, Xe, Xs)
I = diag(1,n_sample)
Hols = X%*%(solve(t(X)%*%X))%*%t(X)
RH0 = t(I-Hols)%*%(I-Hols)
epsilon = (I-Hols)%*%y
#compute R(H1)
I= diag(1,n_sample)
B = sec$B
Xcc = cbind(rep(1,n_sample), Xc)
H0 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H0)%*%(I-H0)
save(RH1, file="RH1_stationary_rotD50pga.RData")
T0 = (t(y) %*% (RH0-RH1) %*% y) / (t(y) %*% RH1 %*% y)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
print("Permutations")
pb = progress_bar$new(total=n_perm, format = "  computing [:bar] :percent eta: :eta")
pb$tick(0)
for (i in 1:n_perm){
  eps_star = sample(epsilon)
  y_star = H0 %*% y + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
  pb$tick()
}
p_cumulative_stationary = sum(t_stat>as.numeric(T0))/n_perm
p_cumulative_stationary

