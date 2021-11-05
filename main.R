
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
source("midpoint/functions/SEC_only_calibration.R")
source("functions/ESC_only_calibration.R")
source("functions/SEC_only_constant_intercept_calibration.R")
source("functions/SEC_no_intercept_calibration.R")
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
    print(paste0("Bandwidth event: ", bwe_tot[i], " / Bandwidth site: ", bws_tot[j]))
    bandwidth = gcv_mei_only_one(bwe       = bwe_tot[i],
                                 bws       = bws_tot[j],
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
## Joint test for the stationarity of the coefficients --------------------
# H0: all coefficients are constant
# H1: at least one coefficient is non-stationary
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

#save(only_intercept, file="res/only_intercept.RData")
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
for (i in 1:n_perm){
  eps_star = sample(epsilon)
  y_star = Hols %*% y + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
  pb$tick()
}
p_ols = sum(t_stat>as.numeric(T0))/n_perm
p_ols

## One-at-a-time test for the stationarity of coefficients ---------------------------
coef_to_check = "b1" #change it among: {"b1","b2","c1","c2","c3","f1","f2","k"}
regs = list(b1 = b1,
            b2 = b2,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            f1 = f1,
            f2 = f2,
            k  = k)
p_stationarity = stationarity_check(coef_to_check = coef_to_check,
                                    regs          = regs,
                                    y             = y,
                                    bwe           = bwe,
                                    bws           = bws,
                                    utm_ev_sp     = utm_ev_sp,
                                    utm_st_sp     = utm_st_sp)
p_stationarity

## Joint test for the stationarity of coefs found stationary -----------------------------------------------------------------------
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
for (i in 1:n_perm){
  eps_star = sample(epsilon)
  y_star = H0 %*% y + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
  pb$tick()
}
p_cumulative_stationary = sum(t_stat>as.numeric(T0))/n_perm
p_cumulative_stationary

## One-at-a-time test for significance of constant coefficients ------------------------------------------
# H0: a constant coefficient is null
# H1: a constant coefficient is different from zero
# Note that: as R(H1) one may take R(H1) from the previous point
Xc = cbind(b1,b2,f1,f2,c1)
names_Xc = c("b1","b2","f1","f2","c1")
Xe = cbind(c2,c3)
Xs = k

# NT: you may change the three matrices above accordingly to the results of the
# your stationarity checks

# WARNING: I don't know what happens if there is no coefficient varying with
#          site or with event (i.e. either Xe or Xs are null)

coef_to_check = "intercept" #change it among all regs that are found constant

p_significance = significance_check(coef_to_check = coef_to_check,
                                    names_Xc      = names_Xc,
                                    Xc            = Xc,
                                    Xe            = Xe,
                                    Xs            = Xs,
                                    y             = y,
                                    bwe           = bwe,
                                    bws           = bws,
                                    utm_ev_sp     = utm_ev_sp,
                                    utm_st_sp     = utm_st_sp)
p_significance


## GCV comparison --------------------------------------------------------------
# Comparison of the GCV values obtained using SEC and ESC

Xc = cbind(b1,b2,f1,f2,c1)
Xe = cbind(c2,c3)
Xs = k

gcvESC = gcv_mei_only_one(bwe,bws,ESC_only_calibration, Xc, Xe, Xs, y, "c", coordinates(utm_ev_sp),
                          coordinates(utm_st_sp))
print("ESC: ", gcvESC)

gcvSEC = gcv_mei_only_one(bwe,bws,SEC_only_calibration, Xc, Xe, Xs, y, "c", coordinates(utm_ev_sp),
                          coordinates(utm_st_sp))
print("SEC: ", gcvSEC)

## Computation of $R^2_{adj}$ ---------------------------------------------------

Xc = cbind(b1,b2,f1,f2,c1)
Xe = cbind(c2,c3)
Xs = k

#sec = SEC_only_calibration(Xc        = Xc,
#                           Xe        = Xe,
#                           Xs        = Xs,
#                           y         = y,
#                           intercept = "c",
#                           bwe       = bwe,
#                           bws       = bws,
#                           utm_ev_sp = coordinates(utm_ev_sp),
#                           utm_st_sp = coordinates(utm_st_sp))

source("parallel/functions/SEC_only_calibration.R")
SEC_only_calibration(Xc        = Xc,
                     Xe        = Xe,
                     Xs        = Xs,
                     y         = y,
                     intercept = "c",
                     bwe       = bwe,
                     bws       = bws,
                     utm_ev_sp = coordinates(utm_ev_sp),
                     utm_st_sp = coordinates(utm_st_sp),
                     model     = "benchmark")

load("benchmark/large_matrices/only_calibration_Hs.RData")
load("benchmark/large_matrices/only_calibration_He.RData")

#create B
N = length(y)
I = diag(rep(1,N))
B = I - He - Hs + Hs %*% He
save(B, file="benchmark/large_matrices/only_calibration_B.RData")
rm(Hs,He)

Xcc = cbind(rep(1,N), Xc)
H = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
save(H, file="benchmark/large_matrices/only_calibration_H.RData")
rm(B)

epsilon= (I-H)%*%y
delta1 = N-2*tr(H)+tr(t(H)%*%H) #this can be saved as "delta1_pga.RData"
delta2 = tr((t(I-H)%*%(I-H)) %*% (t(I-H)%*%(I-H))) #this can be saved as "delta1_pga.RData"
rss = sum(epsilon^2)
sigma2hat = rss/delta1
tss = sum((H%*%y-mean(y))^2)
sqrt(sigma2hat)
R2 = 1-rss/tss
R2adj = 1-(1-R2)*(N-1)/delta1
R2adj

## Calibration  ------------------------------------------
# Computation of the regression coefficients
source("parallel/functions/SEC_grid_creation.R")

result = SEC_grid_creation(Xc        = Xc,
                           Xe        = Xe,
                           Xs        = Xs,
                           y         = y,
                           intercept = "c",
                           bwe       = bwe,
                           bws       = bws,
                           utm_ev_sp = coordinates(utm_ev_sp),
                           utm_st_sp = coordinates(utm_st_sp),
                           grid      = coords_utm,
                           model     = "benchmark")

save(result, file="benchmark/large_matrices/result.RData")

load("benchmark/large_matrices/result.RData")

beta_const = result$beta_c

beta_k = t(result$beta_s)
beta_k_coord = cbind(coords_utm, t(beta_k))
beta_k_coord = as.data.frame(beta_k_coord)
range(beta_k)

beta_c2 = result$beta_e[1,]
beta_c2_coord = cbind(coords_utm, beta_c2)
beta_c2_coord = as.data.frame(beta_c2_coord)
range(beta_c2)

beta_c3 = result$beta_e[2,]
beta_c3_coord = cbind(coords_utm, beta_c3)
beta_c3_coord = as.data.frame(beta_c3_coord)
range(beta_c3)

## PLOTS ---------------------------------------------------------
# Plots of non-stationary regression coefficients

## k
beta_k_bis = beta_k
beta_k_bis[beta_k_bis<(-1)]=-1
beta_k_bis[beta_k_bis>(0)]=0
beta_k_coord = cbind(coords_utm, t(beta_k))
beta_k_coord = as.data.frame(beta_k_coord)
beta_k_bis_coord = cbind(coords_utm, t(beta_k_bis))
beta_k_bis_coord = as.data.frame(beta_k_bis_coord)
x11()
ggplot() + 
  geom_tile(beta_k_bis_coord, mapping = aes(x=x1, y=x2, fill=beta_k_bis))+
  scale_fill_gradientn(colours = c("darkblue", "dodgerblue1", "cadetblue2", "white"), name = "k"
                       ,limits = c(-1,0), breaks = c(-1, -0.5,0)
  )+
  geom_sf(data = shape_utm_no_lamp, size = 1.6, color = "black", fill = NA )+
  geom_point(data = utm_st, aes(x=longitude_st, y=latitude_st), fill= 'firebrick3',
             size = 1, shape = 21, stroke = 0.7)+
  ggtitle("Benchmark")+
  theme(axis.text=element_text(size=15, colour = "black"),
        axis.title=element_blank(),
        legend.text=element_text(size=15, colour = "black"),
        plot.title = element_text(size=15, hjust = 0.5),
        axis.ticks=element_blank(),
        legend.title = element_text(size=15, colour = "black"),
        panel.background = element_rect(fill = "lightcyan2", colour = "skyblue3",
                                        size = 2, linetype = "solid"),
        panel.grid = element_line(size = 0.25, linetype = 'solid',
                                  colour = "aliceblue"))
ggsave(filename = "k.png",
       plot = last_plot(),
       device = NULL,
       path = "benchmark/coefs_estimates",
       scale = 1,
       limitsize = FALSE,
       dpi = 320)

dev.off()

## c2
beta_c2_bis = beta_c2
beta_c2_bis[beta_c2_bis<(-2)]=-2
beta_c2_bis[beta_c2_bis>(-1)]=-1
beta_c2_coord = cbind(coords_utm, beta_c2)
beta_c2_coord = as.data.frame(beta_c2_coord)
beta_c2_bis_coord = cbind(coords_utm, beta_c2_bis)
beta_c2_bis_coord = as.data.frame(beta_c2_bis_coord)
x11()
ggplot() + 
  geom_tile(beta_c2_bis_coord, mapping = aes(x=x1, y=x2, fill=beta_c2_bis))+
  #coord_sf()+
  scale_fill_gradientn(colours = c("darkblue", "dodgerblue1", "cadetblue2", "white")
                       , name = expression(paste(c[2]))
                       ,limits = c(-2,-1), breaks = c(-2, -1.5, -1),
                       labels = c(-2, -1.5, -1))+
  geom_sf(data = shape_utm_no_lamp, size = 1.6, color = "black", fill = NA)+
  geom_point(data = utm_ev, aes(x=longitude_ev, y=latitude_ev), fill= 'firebrick3',
             size = 3, shape = 21, stroke = 1.5)+
  ggtitle("Benchmark")+
  theme(axis.text=element_text(size=15, colour = "black"),
        axis.title=element_blank(),
        legend.text=element_text(size=15, colour = "black"),
        plot.title = element_text(size=15, hjust = 0.5),
        axis.ticks=element_blank(),
        legend.title = element_text(size=15, colour = "black"),
        panel.background = element_rect(fill = "lightcyan2", colour = "skyblue3",
                                        size = 2, linetype = "solid"),
        panel.grid = element_line(size = 0.25, linetype = 'solid',
                                  colour = "aliceblue"))
ggsave(filename = "c2.png",
       plot = last_plot(),
       device = NULL,
       path = "benchmark/coefs_estimates",
       scale = 1,
       limitsize = FALSE,
       dpi = 320)

dev.off()

## c3
beta_c3_bis = beta_c3
beta_c3_bis[beta_c3_bis<(-0.009)]=-0.009
beta_c3_bis[beta_c3_bis>(0.005)]=0.005
beta_c3_coord = cbind(coords_utm, beta_c3)
beta_c3_coord = as.data.frame(beta_c3_coord)
beta_c3_bis_coord = cbind(coords_utm, beta_c3_bis)
beta_c3_bis_coord = as.data.frame(beta_c3_bis_coord)
x11()

ggplot() + 
  geom_tile(beta_c3_bis_coord, mapping = aes(x=x1, y=x2, fill=beta_c3_bis))+
  scale_fill_gradientn(colours = c("darkblue", "dodgerblue1", "cadetblue2", "white"),
                       name = expression(paste(c[3])),
                       limits = c(-0.01, 0.005), breaks = c(-0.01, -0.005, 0, 0.005),
                       labels = c(-0.01, -0.005, 0, 0.005)) +
  geom_sf(data = shape_utm_no_lamp, size = 1.6, color = "black", fill = NA)+
  geom_point(data = utm_ev, aes(x=longitude_ev, y=latitude_ev), fill= 'firebrick3',
             size = 3, shape = 21, stroke = 1.5)+
  ggtitle("Benchmark") +
  theme(axis.text=element_text(size=15, colour = "black"),
        axis.title=element_blank(),
        legend.text=element_text(size=15, colour = "black"),
        plot.title = element_text(size=15, hjust = 0.5),
        axis.ticks=element_blank(),
        legend.title = element_text(size=15, colour = "black"),
        panel.background = element_rect(fill = "lightcyan2", colour = "skyblue3",
                                        size = 2, linetype = "solid"),
        panel.grid = element_line(size = 0.25, linetype = 'solid',
                                  colour = "aliceblue"))

ggsave(filename = "c3.png",
       plot = last_plot(),
       device = NULL,
       path = "benchmark/coefs_estimates",
       scale = 1,
       limitsize = FALSE,
       dpi = 320)

dev.off()