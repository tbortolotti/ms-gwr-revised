#'
#' ESC algorithm for full calibration and computation of the regression coefficients
#'
#' This function builds matrices He, Hs and B that are used for the MS-GWR
#' via ESC algorithm and estimates the regression coefficients over a grid of coordinates given as input
#'
#' @param Xc:          matrix of constant predictor variables
#' @param Xe:          matrix of event-varying predictor variables
#' @param Xs:          matrix of site-varying predictor variables
#' @param y:           response variable
#' @param intercept:   either "c" (constant), "e" (event-dependent), "s" (site-dependent)
#' @param bw1:         bandwidth for event or midpoint
#' @param bw2:         bandwidth for site
#' @param utm_1_sp:    utm coordinates of the events or midpoints
#' @param utm_2_sp:    utm coordinates of the site
#' @param grid:        coordinates of the grid over which to estimate the regression coefficients
#' 
#' @return a six-element list with the following components:
#'         beta_c:      estimate of the constant coefficients
#'         beta_s:      estimate of site-dependent coefficients
#'         beta_e:      estimate of event-dependent coefficients 
#'         He:          matrix He
#'         Hs:          matrix Hs
#'         B:           matrix B
#'

ESC_general = function(Xc, Xe, Xs, y,intercept, bw1, bw2, utm_1_sp, utm_2_sp, grid){
  dist_e_sim = gw.dist(utm_1_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim = gw.dist(utm_2_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_e_sim_cal = gw.dist(utm_1_sp, utm_1_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_2_sp, utm_2_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  beta_c = rep(0,n_c)
  L = dim(grid)[1]
  beta_e = matrix(0,n_e, L)
  beta_s = matrix(0,n_s, L)
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bw1))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    He[i,] = Xe[i,] %*% Ae
  }
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bw2))
    As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
    Hs[i,] = Xs[i,] %*% As
  }
  
  #create B
  B = I - He - Hs + He %*% Hs
  
  #1)compute beta_c
  beta_c = solve(t(Xc)%*%t(B)%*%B%*%Xc)%*%t(Xc)%*%t(B)%*%B%*%y
  
  #2)compute y_tilde
  y_tilde = y-Xc%*%beta_c
  
  #3)compute beta_s for whole grid
  for (i in 1:L){
    Ws = diag(gauss_kernel(dist_s_sim[,i],bw2))
    As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
    beta_s[,i] = As %*% y_tilde
  }
  
  #4)compute y_tilde_s
  y_tilde_s = (I-Hs)%*%y_tilde
  
  #5)compute beta_e for whole grid
  for (i in 1:L){
    We = diag(gauss_kernel(dist_e_sim[,i],bw1))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    beta_e[,i] = Ae %*% y_tilde_s
  }
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e,
                "He" = He,
                "Hs" = Hs,
                "B" = B)
  return(betas)
}