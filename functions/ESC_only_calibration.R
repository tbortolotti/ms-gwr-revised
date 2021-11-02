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
#' 
#' @return a three-element list with the following components:
#'         He:    matrix He
#'         Hs:    matrix Hs
#'         B:     matrix B
#'

ESC_only_calibration = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp){
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
  
  #create He
  print("Create He")
  pb = progress_bar$new(total=N, format = "  computing [:bar] :percent eta: :eta")
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
    He[i,] = Xe[i,] %*% Ae
    #print(c("He",i))
    pb$tick()
  }
  
  #create Hs
  print("Create Hs")
  pb = progress_bar$new(total=N, format = "  computing [:bar] :percent eta: :eta")
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
    Hs[i,] = Xs[i,] %*% As
    #print(c("Hs",i))
    pb$tick()
  }
  
  #create B
  B = I - He - Hs + He %*% Hs
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  
  return(calibration)
}
