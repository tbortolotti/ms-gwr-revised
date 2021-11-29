#'
#' SEC algorithm for calibration
#'
#' This function builds matrices He, Hs and B that are used for the MS-GWR
#' via SEC algorithm
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
#' 
#' @return a three-element list with the following components:
#'         He:    matrix He
#'         Hs:    matrix Hs
#'         B:     matrix B
#'

SEC_calibration = function(Xc, Xe, Xs, y,intercept, bw1, bw2, utm_1_sp, utm_2_sp){
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
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create Hs
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bw2))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    Hs[i,] = Xs[i,] %*% As
  }
  
  #create He
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bw1))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    He[i,] = Xe[i,] %*% Ae
  }
  
  #create B
  B = I - He - Hs + Hs %*% He
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  
  return(calibration)
}
