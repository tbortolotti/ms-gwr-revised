#'
#' Calibration of a model with constant intercept, via SEC algorithm
#'
#' This function performs the calibration of a model with constant intercept,
#' using the SEC algorithm. In the practice, it builds the matrices
#' He, Hs and B that are used to obtain the estimates of the regression coefficients.
#'
#' @param Xc:          matrix of constant predictor variables
#' @param Xe:          matrix of event-varying predictor variables
#' @param Xs:          matrix of site-varying predictor variables
#' @param y:           response variable
#' @param bwe:         bandwidth for event
#' @param bws:         bandwidth for site
#' @param utm_ev_sp:    utm coordinates of the events
#' @param utm_st_sp:   utm coordinates of the site   
#' 
#' @return a three-element list with the following components:
#'         He:    matrix He
#'         Hs:    matrix Hs
#'         B:     matrix B
#' 

SEC_only_constant_intercept_calibration = function(Xe, Xs, y, bwe, bws, utm_ev_sp, utm_st_sp){
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  Xc = rep(1,N)
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = 1
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  Ae = matrix(0,n_e,N)
  As = matrix(0,n_s,N)
  
  He = matrix(0,N,N)
  Hs = matrix(0,N,N)
  
  #create Hs
  print("Create Hs")
  pb = progress_bar$new(total=N, format = "  computing [:bar] :percent eta: :eta")
  pb$tick(0)
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    Hs[i,] = Xs[i,] %*% As
    #print(c("Hs",i))
    pb$tick()
  }
  
  #create He
  print("Create He")
  pb = progress_bar$new(total=N, format = "  computing [:bar] :percent eta: :eta")
  pb$tick(0)
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    He[i,] = Xe[i,] %*% Ae
    #print(c("He",i))
    pb$tick()
  }
  
  #create B
  B = I - He - Hs + Hs %*% He
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  return(calibration)
}
