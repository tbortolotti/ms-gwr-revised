#'
#' Mixed GWR that performs SC algorithm for calibration
#'
#' This function builds matrix Hs to perform a classical Mixed GWR,
#' in which there is only one type of spatial dependence of the regression coefficients
#' estimates
#'
#' @param Xc:          matrix of constant predictor variables
#' @param Xs:          matrix of spatially-varying predictor variables
#' @param y:           response variable
#' @param bw:          bandwidth
#' @param intercept:   either "c" (constant), "s" (non-stationary) or "null" (the model to be calibrated has no intercept)
#' @param utm_sp:      utm coordinates of the points
#' 
#' @return Hs:    matrix Hs of the MGWR
#'

mixed_SC_calibration = function(Xc, Xs, y, bw, intercept, utm_sp){
  
  N = length(y) #y vector of responses
  
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  }
  else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  I = diag(rep(1,N))
  n_c = dim(Xc)[2] #number constant of covariates
  if (is.null(n_c)) {n_c = 1}
  n_s = dim(Xs)[2] #number of site-dependent covariates
  if (is.null(n_s)) {n_s = 1}
  Hs = matrix(0,N,N)
  dist_s_sim_cal = gw.dist(utm_sp, utm_sp, focus=0, p=2, theta=0, longlat=F)
  
  #create Hs
  if (n_s > 1){
    for (i in 1:N){
      Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bw))
      As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
      Hs[i,] = Xs[i,] %*% As
    }
  }
  
  else if (n_s == 1){
    for (i in 1:N){
      Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bw))
      As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
      Hs[i,] = Xs[i] %*% As
    }
  }
  
  betas <- list("Hs" = Hs)
  return(betas)
}