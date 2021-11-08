#'
#' ESC algorithm for calibration
#'
#' This function builds matrices He, Hs and B that are used for the MS-GWR
#' via ESC algorithm
#'
#' @param Xc:          matrix of constant predictor variables
#' @param Xe:          matrix of event-varying predictor variables
#' @param Xs:          matrix of site-varying predictor variables
#' @param y:           response variable
#' @param intercept:   either "c" (constant), "e" (event-dependent), "s" (site-dependent)
#' @param bwe:         bandwidth for event
#' @param bws:         bandwidth for site
#' @param utm_ev_sp:   utm coordinates of the events
#' @param utm_st_sp:   utm coordinates of the site
#' @param model:       choose among ("midpoint","benchmark") or whichever other model you're working with
#' @param test:        It is a name we give to the calibration and it is a parameter that we use to save
#'                     the large matrices generated
#'                     
#' @return a three-element list with the following components:
#'         He:    matrix He
#'         Hs:    matrix Hs
#'         B:     matrix B
#'

ESC_calibration = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, model, test){
  dist_e_sim_cal = gw.dist(utm_ev_sp, utm_ev_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_st_sp, utm_st_sp, focus=0, p=2, theta=0, longlat=F)
  
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
  
  I = diag(rep(1,N))
  
  ## FUNCTIONS ----------------------------------
  gauss_kernel = function(d, h){
    wgts = exp(-0.5*(d/h)^2)
    return(wgts)
  }
  
  my_fun_e = function(i, Xe, dist_e, bwe, gauss_kernel)
  {
    We = diag(gauss_kernel(dist_e[,i],bwe))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    return(Xe[i,] %*% Ae)
  }
  
  my_fun_s = function(i, Xs, He, dist_s, bws, gauss_kernel)
  {
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
    return(Xs[i,] %*% As)
  }
  ## ---------------------------------------------
  
  #create He
  ncpu = 4 # init cluster parallelization
  sfInit(par=TRUE,cp=ncpu)
  reps = 1:N
  (Start.Time <- Sys.time())
  He = sfSapply(x = reps,
                fun = my_fun_e,
                Xe = Xe,
                dist_e = dist_e_sim_cal,
                bwe = bwe,
                gauss_kernel = gauss_kernel)
  End.Time <- Sys.time()
  print(paste0("Building He: ",round(End.Time - Start.Time, 2)))
  sfStop() #stop cluster parallelization
  He = t(He)
  save(He, file=paste0(model,"/large_matrices/ESC_",test,"_calibration_He.RData"))
  
  #create Hs
  sfInit(par=TRUE,cp=ncpu)
  reps = 1:N
  (Start.Time <- Sys.time())
  Hs = sfSapply(x = reps,
                fun = my_fun_s,
                Xs = Xs,
                He = He,
                dist_s = dist_s_sim_cal,
                bws = bws,
                gauss_kernel = gauss_kernel)
  End.Time <- Sys.time()
  print(paste0("Building Hs: ",round(End.Time - Start.Time, 2)))
  sfStop() #stop cluster parallelization
  Hs = t(Hs)
  save(Hs, file=paste0(model,"/large_matrices/ESC_",test,"_calibration_Hs.RData"))
  
  #create B
  B = I - He - Hs + He %*% Hs
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  
  return(calibration)
}
