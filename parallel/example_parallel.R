setwd("C:/Users/Teresa Bortolotti/Documents/R/ms-gwr-revised")

rm(list=ls())
graphics.off()
cat("\014")

# UTILITIES ---------------------------------------------------------
load("parallel/support/utilities_for_parallel_trials.RData")

# FUNCTION ----------------------------------------------------------

gauss_kernel = function(d, h){
  wgts = exp(-0.5*(d/h)^2)
  return(wgts)
}

my_fun = function(i, Xs, dist_mat, bws, gauss_kernel)
{
  Ws = diag(gauss_kernel(dist_mat[,i],bws))
  As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
  return(Xs[i,] %*% As)
}


# LOOP -------------------------------------------------------------------------

library(snowfall)

# init cluster parallelization
ncpu = 2
sfInit(par=TRUE,cp=ncpu)

# number of simulations
N = dim(Xs)[1]
#N = 10
reps = 1:N

(Start.Time <- Sys.time())
Hs = sfSapply(x = reps,
             fun = my_fun,
             Xs = Xs,
             dist_mat = dist_s_sim_cal,
             bws = bws,
             gauss_kernel = gauss_kernel)
End.Time <- Sys.time()
print(paste0("Parallel: ",round(End.Time - Start.Time, 2)))
# stop cluster parallelization
sfStop()
Hs = t(Hs)
save(Hs, file="parallel/support/Hs.RData")

NN = 10
N = dim(Xs)[1]
Hs_vero = matrix(0,N,N)
(Start.Time <- Sys.time())
for (i in 1:N){
  Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
  As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
  Hs_vero[i,] = Xs[i,] %*% As
}
End.Time <- Sys.time()
print(paste0("Serie: ",round(End.Time - Start.Time, 2)))
save(Hs_vero, file="parallel/support/Hs_vero.RData")


sum(Hs!=Hs_vero)

# AIUTAMI -------------------------------------
i=1
Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
aiut <- diag(Ws)

x11(width=8000, height=5000)
hist(log10(aiut))

length(which(aiut<=1e-5))
log10(0.01)

