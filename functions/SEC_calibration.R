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
#' @param model:       choose among ("midpoint","benchmark") or whichever other model you're working with. It is for saving purposes.
#' @param test:        It is a name we give to the calibration and it is a parameter that we use to save
#'                     the large matrices generated
#' 
#' @return a three-element list with the following components:
#'         He:    matrix He
#'         Hs:    matrix Hs
#'         B:     matrix B
#'

SEC_calibration = function(Xc, Xe, Xs, y, intercept, bw1, bw2, utm_1_sp, utm_2_sp, model, test="NULL"){
  
  dist_e_sim_cal = gw.dist(utm_1_sp, utm_1_sp, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim_cal = gw.dist(utm_2_sp, utm_2_sp, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  } else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  } else if (intercept == "s"){
    Xs = cbind(rep(1,N), Xs)
  }
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  I = diag(rep(1,N))
  
  ## FUNCTIONS ----------------------------------
  gauss_kernel = function(d, h){
    wgts = exp(-0.5*(d/h)^2)
    return(wgts)
  }
  
  my_fun_s = function(i, Xs, dist_s, bw2, gauss_kernel)
  {
    Ws = diag(gauss_kernel(dist_s[,i],bw2))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    return(Xs[i,] %*% As)
  }
  
  my_fun_e = function(i, Xe, Hs, dist_e, bw1, gauss_kernel)
  {
    We = diag(gauss_kernel(dist_e[,i],bw1))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    return(Xe[i,] %*% Ae)
  }
  
  ## ---------------------------------------------
  
  #create Hs
  ncpu = 4
  reps = 1:N
  (Start.Time <- Sys.time())
  sfInit(par=TRUE,cp=ncpu)
  Hs = sfSapply(x = reps,
                fun = my_fun_s,
                Xs = Xs,
                dist_s = dist_s_sim_cal,
                bw2 = bw2,
                gauss_kernel = gauss_kernel)
  sfStop()
  End.Time <- Sys.time()
  print(paste0("Building Hs: ",round(End.Time - Start.Time, 2)))
  
  Hs = t(Hs)
  if(test=="full_computation")
  {
    save(Hs, file=paste0(model,"/large_matrices/SEC_",test,"_calibration_Hs.RData"))
  }
  
  #create He
  (Start.Time <- Sys.time())
  sfInit(par=TRUE,cp=ncpu)
  He = sfSapply(x = reps,
                fun = my_fun_e,
                Xe = Xe,
                Hs = Hs,
                dist_e = dist_e_sim_cal,
                bw1 = bw1,
                gauss_kernel = gauss_kernel)
  sfStop()
  End.Time <- Sys.time()
  print(paste0("Building He: ",round(End.Time - Start.Time, 2)))
  He = t(He)
  if(test=="full_computation")
  {
    save(He, file=paste0(model,"/large_matrices/SEC_",test,"_calibration_He.RData"))
  }
  
  #create B
  B = I - He - Hs + Hs %*% He
  if(test=="full_computation")
  {
    save(B, file=paste0(model,"/large_matrices/SEC_",test,"_calibration_B.RData"))
  }
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  
  return(calibration)
  
}

