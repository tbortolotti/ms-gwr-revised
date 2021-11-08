# Working directory (change)
setwd("/home/giovanni/Scrivania/teresa/ms-gwr-revised")

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
library(ggplot2)
library(snowfall)

rm(list=ls())
graphics.off()
cat("\014")

## Load -------------------------------------------------
# Load coordinates data and regressors (SEE file data_preparation.R)
load("data_dir/utm_coordinates.RData")
load("data_dir/regressors.RData")

# Load functions
source("functions/gcv_mei_only_one.R")
source("parallel/functions/ESC_calibration.R")
source("parallel/functions/SEC_calibration.R")
source("functions/stationarity_check.R")
source("functions/significance_check.R")

## Selection of the optimal bandwidth -------------------------------
bwe_tot = c(10000, 25000, 50000, 75000, 100000)
bws_tot = c(10000, 25000, 50000, 75000, 100000)
bwm_tot = c(10000, 25000, 50000, 75000, 100000)
n_bwe = length(bwe_tot)
n_bws = length(bws_tot)
n_bwm = length(bwm_tot)
bw_table = matrix(0, n_bwm, n_bws)
row.names(bw_table) = bwm_tot 
colnames(bw_table) = bws_tot

for (i in 1:n_bwm){ #i=1
  for (j in 1:n_bws){ # j=1
    print(paste0("Bandwidth midpoint: ", bwm_tot[i], " / Bandwidth site: ", bws_tot[j]))
    bandwidth = gcv_mei_only_one(bwe       = bwm_tot[i],
                                 bws       = bws_tot[j],
                                 func      = SEC_calibration,
                                 Xc        = Xc,
                                 Xe        = Xe,
                                 Xs        = Xs,
                                 y         = y,
                                 intercept = "c",
                                 utm_ev_sp = coordinates(utm_md_sp),
                                 utm_st_sp = coordinates(utm_st_sp),
                                 model     = "midpoint")
    bw_table[i,j] = bandwidth$gcv
  }
}



bw_best = which(bw_table == max(bw_table), arr.ind = T)
bwm = bwm_tot[bw_best[1]]
bws = bws_tot[bw_best[2]]

bwm = 25000
bws = 75000

## PERMUTATION TESTS -----------------------------------------------------------
## Joint test for the stationarity of the coefficients --------------------
# H0: all coefficients are constant
# H1: at least one coefficient is non-stationary
Xc = cbind(b1,b2,f1,f2,c1)
Xe = cbind(c2,c3)
Xs = k

ols = lm(y ~ Xc + Xe + Xs)

# only the intercept is considered as constant, hence Xc is given to SEC_calibration function as empty
Xc = c()
Xe = cbind(c2,c3)
Xs = cbind(b1,b2,f1,f2,c1,k)

(Start.Time <- Sys.time())
only_intercept = SEC_calibration(Xc        = Xc,
                                 Xe        = Xe,
                                 Xs        = Xs,
                                 intercept = "c",
                                 y         = y,
                                 bwe       = bwm,
                                 bws       = bws,
                                 utm_ev_sp = coordinates(utm_md_sp),
                                 utm_st_sp = coordinates(utm_st_sp),
                                 model     = "midpoint",
                                 test      = "only_constant_intercept")

End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

N = length(y)
#compute R(H0)
X = cbind(rep(1,N), Xc, Xe, Xs)
I = diag(1,N)
Hols = X%*%(solve(t(X)%*%X))%*%t(X)
RH0 = t(I-Hols)%*%(I-Hols)
epsilon = (I-Hols)%*%y
#compute R(H1)
B = only_intercept$B
Xcc = rep(1,N)
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
save(RH1, file="midpoint/RH1_only_intercept_rotD50pga.RData")
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
p = sum(t_stat>as.numeric(T0))/n_perm
p
save(p, file="midpoint/pvals/model_stationarity.RData")

## One-at-a-time test for the stationarity of coefficients ---------------------------
n_coef_to_check = "b2" #change it among: {"b1","b2","c1","c2","c3","f1","f2","k"}
coef_to_check = b2
regs = list(b1 = b1,
            b2 = b2,
            c1 = c1,
            c2 = c2,
            c3 = c3,
            f1 = f1,
            f2 = f2,
            k  = k)
(Start.Time <- Sys.time())
p = stationarity_check(n_coef_to_check = n_coef_to_check,
                       coef_to_check   = coef_to_check,
                       regs            = regs,
                       y               = y,
                       bwe             = bwm,
                       bws             = bws,
                       utm_ev_sp       = utm_md_sp,
                       utm_st_sp       = utm_st_sp,
                       model           = "midpoint")
p
save(p, file=paste0("midpoint/pvals/stationarity_",n_coef_to_check,".RData"))
End.Time <- Sys.time()
round(End.Time - Start.Time, 2)

## Joint test for the stationarity of coefs found stationary -----------------------------------------------------------------------
# H0: b1, b2, f1, f2, c1, other than the intercept, are constant
# H1: at least one is spatially varying
Xc = cbind(b1,b2,f1,f2,c1)
Xe = cbind(c2,c3)
Xs = k
ols = lm(y ~ Xc + Xe + Xs)
sec = SEC_calibration(Xc, Xe, Xs, y, "c", bwe, bws, coordinates(utm_ev_sp), coordinates(utm_st_sp))
#compute R(H0)
X = cbind(rep(1,N), Xc, Xe, Xs)
I = diag(1,N)
Hols = X%*%(solve(t(X)%*%X))%*%t(X)
RH0 = t(I-Hols)%*%(I-Hols)
epsilon = (I-Hols)%*%y
#compute R(H1)
I= diag(1,N)
B = sec$B
Xcc = cbind(rep(1,N), Xc)
H0 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H0)%*%(I-H0)
save(RH1, file="midpoint/RH1_stationary_rotD50pga.RData")
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
p = sum(t_stat>as.numeric(T0))/n_perm
p
save(p, file="midpoint/pvals/joint_stationarity.RData")

## One-at-a-time test for significance of constant coefficients ------------------------------------------
# H0: a constant coefficient is null
# H1: a constant coefficient is different from zero
# Note that: as R(H1) one may take R(H1) from the previous point
Xc = cbind(b1,b2,f1,f2,c1)
names_Xc = c("b1","b2","f1","f2","c1")
Xe = cbind(c2,c3)
Xs = k

# NT: you may change the three matrices above accordingly to the results of
# your stationarity checks

# WARNING: If your stationarity tests result in no site-varying coefficients or no event-varying
# coefficients, than you resort to a single-effect GWR. You can't use the code below, but you
# have to rely o a single-effect GWR. NOTA PER ME: capire come funziona una regressione con
# un solo coefficiente variabile e sviluppare la fuzione accordingly

coef_to_check = "intercept" #change it among all regs that are found constant

p = significance_check(coef_to_check = coef_to_check,
                       names_Xc      = names_Xc,
                       Xc            = Xc,
                       Xe            = Xe,
                       Xs            = Xs,
                       y             = y,
                       bwe           = bwm,
                       bws           = bws,
                       utm_ev_sp     = utm_md_sp,
                       utm_st_sp     = utm_st_sp,
                       model         = "midpoint")
p
save(p, file=paste0("midpoint/pvals/significance_",coef_to_check,".RData"))


## GCV comparison --------------------------------------------------------------
# Comparison of the GCV values obtained using SEC and ESC

Xc = cbind(b1,b2,f1,f2,c1)
Xe = cbind(c2,c3)
Xs = k

gcvESC = gcv_mei_only_one(bwm,bws,ESC_calibration, Xc, Xe, Xs, y, "c", coordinates(utm_md_sp),
                          coordinates(utm_st_sp))
print("ESC: ", gcvESC)

gcvSEC = gcv_mei_only_one(bwm,bws,SEC_calibration, Xc, Xe, Xs, y, "c", coordinates(utm_md_sp),
                          coordinates(utm_st_sp))
print("SEC: ", gcvSEC)

## Computation of $R^2_{adj}$ ---------------------------------------------------

Xc = cbind(b1,b2,f1,f2,c1)
Xe = cbind(c2,c3)
Xs = k

sec = SEC_calibration(Xc        = Xc,
                      Xe        = Xe,
                      Xs        = Xs,
                      y         = y,
                      intercept = "c",
                      bwe       = bwm,
                      bws       = bws,
                      utm_ev_sp = coordinates(utm_md_sp),
                      utm_st_sp = coordinates(utm_st_sp),
                      model     = "midpoint",
                      test      = "full_computation")

load("midpoint/large_matrices/full_calibration_Hs.RData")
load("midpoint/large_matrices/full_calibration_He.RData")

#create B
N = length(y)
I = diag(rep(1,N))
B = I - He - Hs + Hs %*% He
save(B, file="midpoint/large_matrices/full_calibration_B.RData")
rm(Hs,He)

Xcc = cbind(rep(1,N), Xc)
H = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
save(H, file="midpoint/large_matrices/full_calibration_H.RData")
rm(B)
epsilon= (I-H)%*%y
delta1 = N-2*tr(H)+tr(t(H)%*%H)
delta2 = tr((t(I-H)%*%(I-H)) %*% (t(I-H)%*%(I-H)))
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
                           bwe       = bwm,
                           bws       = bws,
                           utm_ev_sp = coordinates(utm_md_sp),
                           utm_st_sp = coordinates(utm_st_sp),
                           grid      = coords_utm,
                           model     = "midpoint")

save(result, file="midpoint/large_matrices/SEC_grid_creation.RData")
load("midpoint/large_matrices/SEC_grid_creation.RData")

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
                       ,limits = c(-1,0), breaks = c(-1, -0.5,0))+
  geom_sf(data = shape_utm_no_lamp, size = 1.6, color = "black", fill = NA )+
  geom_point(data = utm_st, aes(x=longitude_st, y=latitude_st), fill= 'firebrick3',
             size = 1, shape = 21, stroke = 0.7)+
  ggtitle("Midpoint") +
  theme(axis.text=element_text(size=15, colour = "black"),
        axis.title=element_blank(),
        legend.text=element_text(size=15, colour = "black"),
        plot.title = element_text(size=15, hjust=0.5),
        axis.ticks=element_blank(),
        legend.title = element_text(size=15, colour = "black"),
        panel.background = element_rect(fill = "lightcyan2", colour = "skyblue3",
                                        size = 2, linetype = "solid"),
        panel.grid = element_line(size = 0.25, linetype = 'solid',
                                  colour = "aliceblue"))
ggsave(filename = "k.png",
       plot = last_plot(),
       device = NULL,
       path = "midpoint/coefs_estimates",
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
                       ,limits = c(-2, -1), breaks = c(-2, -1.5, -1),
                       labels = c(-2, -1.5, -1))+
  geom_sf(data = shape_utm_no_lamp, size = 1.6, color = "black", fill = NA)+
  geom_point(data = utm_ev, aes(x=longitude_ev, y=latitude_ev), fill= 'firebrick3',
             size = 3, shape = 21, stroke = 1.5)+
  ggtitle("Midpoint") +
  theme(axis.text=element_text(size=15, colour = "black"),
        axis.title=element_blank(),
        legend.text=element_text(size=15, colour = "black"),
        plot.title = element_text(size=15, hjust=0.5),
        axis.ticks=element_blank(),
        legend.title = element_text(size=15, colour = "black"),
        panel.background = element_rect(fill = "lightcyan2", colour = "skyblue3",
                                        size = 2, linetype = "solid"),
        panel.grid = element_line(size = 0.25, linetype = 'solid',
                                  colour = "aliceblue"))
ggsave(filename = "c2.png",
       plot = last_plot(),
       device = NULL,
       path = "midpoint/coefs_estimates",
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
  scale_fill_gradientn(colours = c("darkblue", "dodgerblue1", "cadetblue2", "white")
                       , name = expression(paste(c[3]))
                       , limits = c(-0.01, 0.005), breaks = c(-0.01, -0.005, 0, 0.005),
                       labels = c(-0.01, -0.005, 0, 0.005)) +
  geom_sf(data = shape_utm_no_lamp, size = 1.6, color = "black", fill = NA)+
  geom_point(data = utm_ev, aes(x=longitude_ev, y=latitude_ev), fill= 'firebrick3',
             size = 3, shape = 21, stroke = 1.5)+
  ggtitle("Midpoint")+
  theme(axis.text=element_text(size=15, colour = "black"),
        axis.title=element_blank(),
        legend.text=element_text(size=15, colour = "black"),
        plot.title = element_text(size=15, hjust=0.5),
        axis.ticks=element_blank(),
        legend.title = element_text(size=15, colour = "black"),
        panel.background = element_rect(fill = "lightcyan2", colour = "skyblue3",
                                        size = 2, linetype = "solid"),
        panel.grid = element_line(size = 0.25, linetype = 'solid',
                                  colour = "aliceblue"))

ggsave(filename = "c3.png",
       plot = last_plot(),
       device = NULL,
       path = "midpoint/coefs_estimates",
       scale = 1,
       limitsize = FALSE,
       dpi = 320)

dev.off()