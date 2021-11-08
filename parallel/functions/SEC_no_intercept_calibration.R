#'
#' Calibration of a model without intercept, via SEC algorithm
#'
#' This function performs the calibration of a model without intercept,
#' using the SEC algorithm. In the practice, it builds the matrices
#' He, Hs and B that are used to obtain the estimates of the regression coefficients.
#'
#' @param Xc:          matrix of constant predictor variables
#' @param Xe:          matrix of event-varying predictor variables
#' @param Xs:          matrix of site-varying predictor variables
#' @param y:           response variable
#' @param bwe:         bandwidth for event
#' @param bws:         bandwidth for site
#' @param utm_ev_sp:   utm coordinates of the events
#' @param utm_st_sp:   utm coordinates of the site   
#' @param model:       choose among ("midpoint","benchmark") or whichever other model you're working with
#' 
#' @return a three-element list with the following components:
#'         He:    matrix He
#'         Hs:    matrix Hs
#'         B:     matrix B
#'

SEC_no_intercept_calibration = function(Xc, Xe, Xs, y, bwe, bws, utm_ev_sp, utm_st_sp, model){
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  
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
  
  ## FUNCTIONS ----------------------------------
  gauss_kernel = function(d, h){
    wgts = exp(-0.5*(d/h)^2)
    return(wgts)
  }
  
  my_fun_s = function(i, Xs, dist_s, bws, gauss_kernel)
  {
    Ws = diag(gauss_kernel(dist_s[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    return(Xs[i,] %*% As)
  }
  
  my_fun_e = function(i, Xe, Hs, dist_e, bwe, gauss_kernel)
  {
    We = diag(gauss_kernel(dist_e[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    return(Xe[i,] %*% Ae)
  }
  
  ## ---------------------------------------------
  
  #create Hs
  ncpu = 4 # init cluster parallelization
  sfInit(par=TRUE,cp=ncpu)
  reps = 1:N
  (Start.Time <- Sys.time())
  Hs = sfSapply(x = reps,
                fun = my_fun_s,
                Xs = Xs,
                dist_s = dist_s_sim_cal,
                bws = bws,
                gauss_kernel = gauss_kernel)
  End.Time <- Sys.time()
  print(paste0("Building Hs: ",round(End.Time - Start.Time, 2)))
  sfStop() #stop cluster parallelization
  Hs = t(Hs)
  save(Hs, file=paste0(model,"/large_matrices/no_intercept_calibration_Hs.RData"))
  
  #create He
  sfInit(par=TRUE,cp=ncpu)
  reps = 1:N
  (Start.Time <- Sys.time())
  He = sfSapply(x = reps,
                fun = my_fun_e,
                Xe = Xe,
                Hs = Hs,
                dist_e = dist_e_sim_cal,
                bwe = bwe,
                gauss_kernel = gauss_kernel)
  End.Time <- Sys.time()
  print(paste0("Building He: ",round(End.Time - Start.Time, 2)))
  sfStop() #stop cluster parallelization
  He = t(He)
  save(He, file=paste0(model,"/large_matrices/no_intercept_calibration_He.RData"))
  
  #create B
  B = I - He - Hs + Hs %*% He
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  
  return(calibration)
}