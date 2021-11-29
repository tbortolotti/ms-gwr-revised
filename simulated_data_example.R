# Working directory (change)
setwd("/home/giovanni/Scrivania/teresa/ms-gwr-revised")

# source('package_install.R')

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

## Load functions
source("functions/simulation/functions.R")
source("functions/simulation/bw_cv.R")
source("functions/simulation/SEC_calibration.R")
source("functions/simulation/ESC_calibration.R")
source("functions/simulation/emgwr_calibration_prediction.R")
source("functions/simulation/mixed_SC_calibration.R")
source("functions/simulation/ESC_general.R")
source("functions/simulation/SEC_general.R")

## Building the regression coefficients --------------------------------------------------------------
## Two grids in the range -5 5 are built in order to have the coordinates of station and event
inf = -5
sup = 5
step = 0.5
x = s1 = s2 = e1 = e2 = seq(inf, sup, by = step)
adj = (length(x)+1)/2

grid_e_sim = matrix(0,(length(e1)*length(e2)),2)
for (i in 1:length(e1)){
  for (j in 1:length(e2)){
    grid_e_sim[j+(i-1)*length(s2),] = c(e1[i],e2[j])
  }
}

grid_s_sim = matrix(0,(length(s1)*length(s2)),2)
for (i in 1:length(s1)){
  for (j in 1:length(s2)){
    grid_s_sim[j+(i-1)*length(s2),] = c(s1[i],s2[j])
  }
}

#The regression coefficients are set and, in the case of spatial non-stationarity, they are also plotted----
#beta_c intercept
beta_c_intercept_true = 8

#beta_c variable
beta_c_variable_true = 4

png(file=paste0("plots/simulation/true_coefficients.png"), width=8000, height=5000, units="px", res=800)
par(mfrow=c(1,2), mai = c(0.4, 0.4, 0.4, 0.4))

#beta_e variable
beta_e_variable_true = matrix(0,length(e1), length(e2))
beta_e_variable_true_col = matrix(0,length(e1)*length(e2),3)
for (i in 1:length(e1)){
  for (j in 1:length(e2)){
    #beta_e_variable_true[i,j] = -0.01*(s1[i]+0.5)^3 - 0.01*(s2[j]-0.5)^3 + 0.007*(s2[j]-2.5)^3+3.5
    beta_e_variable_true[i,j] = +0.05*(e1[i]+1)^2 + 0.05*(e2[j])^2 + 1.5 + 0.1*(e1[i]-2)
    beta_e_variable_true_col[j+(i-1)*length(e2),] = c(beta_e_variable_true[i,j],e1[i],e2[j])
  }
}

persp3D(e1,e2,beta_e_variable_true, colvar = beta_e_variable_true,theta=20,
        phi=20,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        #contour = list(col = "grey", side = c("zmin", "z")),
        main=expression(paste(beta["1E"]))
)

#beta_s variable
beta_s_variable_true = matrix(0,length(s1), length(s2))
beta_s_variable_true_col = matrix(0,length(s1)*length(s2),3)
for (i in 1:length(s1)){
  for (j in 1:length(s2)){
    #beta_s_variable_true[i,j] = -0.05*(s1[i]-1)^2 - 0.05*(s2[j])^2 + 4.5 - 0.1*(s1[i]+3)
    beta_s_variable_true[i,j] = -0.01*(s1[i]+0.5)^3 - 0.01*(s2[j]-0.5)^3 + 0.007*(s2[j]-2.5)^3+3.5
    beta_s_variable_true_col[j+(i-1)*length(s2),] = c(beta_s_variable_true[i,j],s1[i],s2[j])
  }
}

persp3D(s1,s2,beta_s_variable_true, colvar = beta_s_variable_true,theta=20,
        phi=20,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        #contour = list(col = "grey", side = c("zmin", "z")),
        main=expression(paste(beta["1S"]))
)
dev.off()

## Building covariates ----------------------------------------------------------------------------------------
## Three covariates are built:
## - stationary covariate
## - site-dependent covariate
## - event-dependent covariate
## After building a random error term, the responses of the model are computed

#create X
set.seed(6494)
n = length(s1)^4
Xc_var = rnorm(n, 3, 5)
Xe_var = rnorm(n, -4, 5)
Xs_var = rnorm(n, 4, 5)

Xc = cbind(rep(1,n), Xc_var)
Xe = Xe_var
Xs = Xs_var

#create variance on four dimensional grid
epsilon = rnorm(n,0,2) 

n = length(s1)^4 #create "four-dimensional grid"
f = 4
l = rep(list(x), f)
coords = expand.grid(l)
X_coords = cbind(Xc_var, Xe_var, Xs_var, coords, epsilon)
names(X_coords)[4:7]=c("E1", "E2", "S1", "S2")

#response y
beta_c_true = cbind(beta_c_intercept_true, beta_c_variable_true)
beta_s_true = beta_s_variable_true_col
beta_e_true = beta_e_variable_true_col
yc = rep(0,n)
for (i in 1:n){
  coords_e = X_coords[i,4:5]
  coords_s = X_coords[i,6:7]
  yc[i] = Xc[i,]%*%t(beta_c_true) +
    Xe_var[i]*beta_e_variable_true[coords_e[[1]]/step+adj,coords_e[[2]]/step+adj] +
    Xs_var[i]*beta_s_variable_true[coords_s[[1]]/step+adj,coords_s[[2]]/step+adj] +
    epsilon[i]
}


#The next step is to create our simulation subset, made of 80 elements picked casually from our generated dataset----
set.seed(2803)
N = 80
lines = sample(1:n, N)
ordered_lines = sort(lines)
X_sim = X_coords[ordered_lines,]
y_sim_c = yc[ordered_lines]

grid_sim = matrix(0,(length(s1)*length(s2)),2)
for (i in 1:length(s1)){
  for (j in 1:length(s2)){
    grid_sim[j+(i-1)*length(s2),] = c(s1[i],s2[j])
  }
}

coords_e_sim = cbind(X_sim$E1, X_sim$E2)
coords_s_sim = cbind(X_sim$S1, X_sim$S2)

## Selection of the optimal bandwidth --------------------------------------------------------------
## Notice that here we always check models with bwe=bws
## In this context, we also check which one is best to be used between ESC and SEC algorithms.

bw_sec = bw_cv(bw_min    = 0.5,
               bw_max    = 3.2,
               step      = 0.15,
               f         = 80,
               func      = SEC_calibration,
               method    = "sec",
               Xc        = X_sim$Xc_var,
               Xe        = X_sim$Xe_var,
               Xs        = X_sim$Xs_var,
               y         = y_sim_c,
               intercept = "c",
               utm_1_sp  = coords_e_sim,
               utm_2_sp  = coords_s_sim)

bw_esc = bw_cv(bw_min    = 0.5,
               bw_max    = 3.2,
               step      = 0.15,
               f         = 80,
               func      = ESC_calibration,
               method    = "esc",
               Xc        = X_sim$Xc_var,
               Xe        = X_sim$Xe_var,
               Xs        = X_sim$Xs_var,
               y         = y_sim_c,
               intercept = "c",
               utm_1_sp  = coords_e_sim,
               utm_2_sp  = coords_s_sim)

cvsum_sec = bw_sec[[2]][,2]
cvsum_esc = bw_esc[[2]][,2]
bws = bw_esc[[2]][,1]
dataset = as.data.frame(cbind(bws, cvsum_sec, cvsum_esc))
melted = melt(dataset, id.vars = "bws")

x11()
ggplot(melted, aes(x=bws, y=value, col=variable)) + 
  geom_line(size = 1.7)+
  labs(y = "CVSS", x = "Bandwidth")+
  scale_color_manual(values=c("royalblue3","skyblue"),labels = c("SEC", "ESC"))+
  theme(axis.text=element_text(size=15, colour = "black"),
        axis.title=element_text(size=25),
        legend.text=element_text(size=25, colour = "black"),
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA)
  )+
  scale_x_continuous(breaks = seq(from = 0.65, to = 3.05, by = 0.3)) +
  geom_vline(xintercept = 1.25, size = 1.2, colour="red3")+
  coord_cartesian(xlim = c(0.6, 3.1),ylim = c(800, 3400))+
  geom_point(size = 4)

ggsave(filename = "bandwidth_cv.png",
       plot = last_plot(),
       device = NULL,
       path = "plots/simulation",
       scale = 1,
       limitsize = FALSE,
       dpi = 320)

dev.off()

## Set optimal bandiwdths
bwe = bws = 1.25

## Permutation tests for spatial non-stationarity --------------------------------------------------------------
# OLS vs completely varying model
ols = lm(y_sim_c ~ X_sim$Xc_var + X_sim$Xe_var + X_sim$Xs_var)


esc_only_intercept = ESC_calibration(Xc        = c(),
                                     Xe        = cbind(X_sim$Xc_var,X_sim$Xe_var),
                                     Xs        = X_sim$Xs_var,
                                     y         = y_sim_c,
                                     intercept = "c",
                                     bw1       = bwe,
                                     bw2       = bws,
                                     utm_1_sp  = coords_e_sim,
                                     utm_2_sp  = coords_s_sim)

#compute R(H0)
X = cbind(rep(1,N), X_sim$Xc_var, X_sim$Xe_var, X_sim$Xs_var)
I = diag(1,N)
Hols = X%*%(solve(t(X)%*%X))%*%t(X)
RH0 = t(I-Hols)%*%(I-Hols)
eps = (I-Hols)%*%y_sim_c
#compute R(H1)
B = esc_only_intercept$B
Xcc = rep(1,dim(X_sim)[1])
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)     
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
set.seed(6494)
pb = progress_bar$new(total=n_perm, format = "  computing [:bar] :percent eta: :eta")
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = Hols %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
  pb$tick()
}
p = sum(t_stat>as.numeric(T0))/n_perm
p

# Non-stationarity of Xc
# We indeed expect it to be constant, i.e. to obtain a high p-value

#calibrate models under H0 and H1
esc = ESC_calibration(Xc        = X_sim$Xc_var,
                      Xe        = X_sim$Xe_var,
                      Xs        = X_sim$Xs_var,
                      y         = y_sim_c,
                      intercept = "c",
                      bw1       = bwe,
                      bw2       = bws,
                      utm_1_sp  = coords_e_sim,
                      utm_2_sp  = coords_s_sim)


esc_only_intercept = ESC_calibration(Xc        = c(),
                                     Xe        = cbind(X_sim$Xc_var,X_sim$Xe_var),
                                     Xs        = X_sim$Xs_var,
                                     y         = y_sim_c,
                                     intercept = "c",
                                     bw1       = bwe,
                                     bw2       = bws,
                                     utm_1_sp  = coords_e_sim,
                                     utm_2_sp  = coords_s_sim)

#compute R(H0)
I= diag(1,N)
B = esc$B
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H0 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH0 = t(I-H0)%*%(I-H0)
eps = (I-H0)%*%y_sim_c
#compute R(H1)
B = esc_only_intercept$B
Xcc = rep(1,dim(X_sim)[1])
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
set.seed(6494)
pb = progress_bar$new(total=n_perm, format = "  computing [:bar] :percent eta: :eta")
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = H0 %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
  pb$tick()
}
p_xc = sum(t_stat>as.numeric(T0))/n_perm
p_xc

# Non-stationarity of Xe
#calibrate models under H0 and H1
esc = ESC_calibration(Xc        = X_sim$Xe_var,
                      Xe        = X_sim$Xc_var,
                      Xs        = X_sim$Xs_var,
                      y         = y_sim_c,
                      intercept = "c",
                      bw1       = bwe,
                      bw2       = bws,
                      utm_1_sp  = coords_e_sim,
                      utm_2_sp  = coords_s_sim)


esc_only_intercept = ESC_calibration(Xc        = c(),
                                     Xe        = cbind(X_sim$Xc_var,X_sim$Xe_var),
                                     Xs        = X_sim$Xs_var,
                                     y         = y_sim_c,
                                     intercept = "c",
                                     bw1       = bwe,
                                     bw2       = bws,
                                     utm_1_sp  = coords_e_sim,
                                     utm_2_sp  = coords_s_sim)

#compute R(H0)
I= diag(1,N)
B = esc$B
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xe_var)
H0 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH0 = t(I-H0)%*%(I-H0)
eps = (I-H0)%*%y_sim_c
#compute R(H1)
B = esc_only_intercept$B
Xcc = rep(1,dim(X_sim)[1])
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
pb = progress_bar$new(total=n_perm, format = "  computing [:bar] :percent eta: :eta")
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = H0 %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
  pb$tick()
}
p_xe = sum(t_stat>as.numeric(T0))/n_perm
p_xe

# Non-stationarity of Xs
# calibrate models under H0 and H1
esc = mixed_SC_calibration(Xc           = X_sim$Xs_var,
                           Xs           = cbind(X_sim$Xe_var, X_sim$Xc_var),
                           y            = y_sim_c,
                           bw           = 1.25,
                           intercept    = "c",
                           utm_sp       = coords_e_sim)

esc_only_intercept = ESC_calibration(Xc        = c(),
                                     Xe        = cbind(X_sim$Xe_var, X_sim$Xc_var),
                                     Xs        = X_sim$Xs_var,
                                     y         = y_sim_c,
                                     intercept = "c",
                                     bw1       = bwe,
                                     bw2       = bws,
                                     utm_1_sp  = coords_e_sim,
                                     utm_2_sp  = coords_s_sim)

#compute R(H0)
I= diag(1,N)
Hs = esc$Hs
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xs_var)
H0 = Hs + (I-Hs)%*%Xcc%*%solve(t(Xcc)%*%t(I-Hs)%*%(I-Hs)%*%Xcc)%*%t(Xcc)%*%t(I-Hs)%*%(I-Hs)
RH0 = t(I-H0)%*%(I-H0)
eps = (I-H0)%*%y_sim_c
#compute R(H1)
B = esc_only_intercept$B
Xcc = rep(1,dim(X_sim)[1])
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
pb = progress_bar$new(total=n_perm, format = "  computing [:bar] :percent eta: :eta")
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = H0 %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
  pb$tick()
}
p_xs = sum(t_stat>as.numeric(T0))/n_perm
p_xs

## Permutation tests for the significance of stationary coefficients ---------------------------------------------
# Check if the intercept is significantly different from 0
## Qui sono quasi sicura ci sia un errore nella definizione dei modelli sotto H0 e H1, ma va beh
## me ne preoccupo solo se viene fuori che il modello è giusto
#calibrate models under H0 and H1
esc_null = mixed_SC_calibration(Xc           = X_sim$Xc_var,
                                Xs           = X_sim$Xs_var,
                                y            = y_sim_c,
                                bw           = 1.25,
                                intercept    = "null",
                                utm_sp       = coords_s_sim)

esc = ESC_calibration(Xc        = c(),
                      Xe        = cbind(X_sim$Xe_var, X_sim$Xc_var),
                      Xs        = X_sim$Xs_var,
                      y         = y_sim_c,
                      intercept = "c",
                      bw1       = bwe,
                      bw2       = bws,
                      utm_1_sp  = coords_e_sim,
                      utm_2_sp  = coords_s_sim)

#compute R(H0)
I= diag(1,N)
Hs = esc_null$Hs
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H0 = Hs + (I-Hs)%*%Xcc%*%solve(t(Xcc)%*%t(I-Hs)%*%(I-Hs)%*%Xcc)%*%t(Xcc)%*%t(I-Hs)%*%(I-Hs)
RH0 = t(I-H0)%*%(I-H0)
eps = (I-H0)%*%y_sim_c
#compute R(H1)
I= diag(1,N)
B = esc$B
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
pb = progress_bar$new(total=n_perm, format = "  computing [:bar] :percent eta: :eta")
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = H0 %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
  pb$tick()
}
p_intercept_null = sum(t_stat>as.numeric(T0))/n_perm
p_intercept_null

#check if Xc is null
## Anche questo penso che sia sbagliato ma nel caso anche qui ci penso a posteriori
#calibrate models under H0 and H1

esc_null = mixed_SC_calibration(Xc           = c(),
                                Xs           = X_sim$Xe_var,
                                y            = y_sim_c,
                                bw           = 1.25,
                                intercept    = "c",
                                utm_sp       = coords_e_sim)

esc = ESC_calibration(Xc        = X_sim$Xc_var,
                      Xe        = X_sim$Xe_var,
                      Xs        = X_sim$Xs_var,
                      y         = y_sim_c,
                      intercept = "c",
                      bw1       = bwe,
                      bw2       = bws,
                      utm_1_sp  = coords_e_sim,
                      utm_2_sp  = coords_s_sim)

#compute R(H0)
I= diag(1,N)
Hs = esc_null$Hs
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H0 = Hs + (I-Hs)%*%Xcc%*%solve(t(Xcc)%*%t(I-Hs)%*%(I-Hs)%*%Xcc)%*%t(Xcc)%*%t(I-Hs)%*%(I-Hs)
RH0 = t(I-H0)%*%(I-H0)
eps = (I-H0)%*%y_sim_c
#compute R(H1)
I= diag(1,N)
B = esc$B
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H1 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
RH1 = t(I-H1)%*%(I-H1)
#compute T
T0 = (t(y_sim_c) %*% (RH0-RH1) %*% y_sim_c) / (t(y_sim_c) %*% RH1 %*% y_sim_c)
#permutations
n_perm = 1000
t_stat = rep(0,n_perm)
pb = progress_bar$new(total=n_perm, format = "  computing [:bar] :percent eta: :eta")
for (i in 1:n_perm){
  eps_star = sample(eps)
  y_star = H0 %*% y_sim_c + eps_star
  t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
  pb$tick()
}
p_xc_null = sum(t_stat>as.numeric(T0))/n_perm
p_xc_null

## Full calibration of the model ---------------------------------------------------------------------------------

bwe = bws = 1.25
esc = ESC_general(Xc        = X_sim$Xc_var,
                  Xe        = X_sim$Xe_var,
                  Xs        = X_sim$Xs_var,
                  y         = y_sim_c,
                  intercept = "c",
                  bw1       = bwe,
                  bw2       = bws,
                  utm_1_sp  = coords_e_sim,
                  utm_2_sp  = coords_s_sim,
                  grid      = as.matrix(grid_sim))

betas_c = esc$beta_c
betas_c

#The spatially varying coefficients can be graphically compared with the true ones----
beta_s = esc$beta_s
beta_e = esc$beta_e

beta_s_grid = matrix(0, length(s1),length(s2))
k=1
for (i in 1:length(s1)){
  for (j in 1:length(s2)){
    beta_s_grid[i,j] = beta_s[1,k]
    k = k+1
  }
}

beta_e_grid = matrix(0, length(e1),length(e2))
k=1
for (i in 1:length(e1)){
  for (j in 1:length(e2)){
    beta_e_grid[i,j] = beta_e[1,k]
    k = k+1
  }
}

png(file=paste0("plots/simulation/coefficient_estimates.png"), width=8000, height=5000, units="px", res=800)
par(mfrow=c(2,2),  mai = c(0.4, 0.4, 0.4, 0.4))
persp3D(e1,e2,beta_e_variable_true, colvar = beta_e_variable_true,theta=30,
        phi=40,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        main=expression(paste(beta["1E"]))
)

persp3D(s1,s2,beta_s_variable_true, colvar = beta_s_variable_true,theta=30,
        phi=40,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        main=expression(paste(beta["1S"]))
)

persp3D(e1,e2,beta_e_grid, colvar = beta_e_grid,theta=30,
        phi=40,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        main=expression(paste(beta["1E, EMGWR"]))
)

persp3D(s1,s2,beta_s_grid, colvar = beta_s_grid,theta=30,
        phi=40,axes= TRUE,scale=2, box=TRUE, nticks=5,
        ticktype="detailed",
        xlab="x", 
        ylab="y",
        zlab="z",
        zlim = c(0,5),
        clim = c(0,5),
        cex.main = 2,
        expand = 0.6,
        bty = "g",
        main=expression(paste(beta["1S, EMGWR"]))
)
dev.off()

## Computation of R^2_{adj} ----------------------------------------------------------------------
I= diag(1,N)
B = esc$B
Xcc = cbind(rep(1,dim(X_sim)[1]), X_sim$Xc_var)
H_esc = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
res = (I-H_esc)%*%y_sim_c
delta1 = N-2*tr(H_esc)+tr(t(H_esc)%*%H_esc)
rss = sum(res^2)
sigma2hat = rss/delta1
tss = sum((H_esc%*%y_sim_c-mean(y_sim_c))^2)
R2 = 1-rss/tss
R2adj = 1-(1-R2)*(N-1)/delta1
R2adj