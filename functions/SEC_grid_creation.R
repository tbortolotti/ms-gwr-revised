#'
#' SEC algorithm for the computation of regression coefficients.
#'
#' This function builds the regression coefficients, each of which evaluated on a
#' grid of points of interest, according to the SEC algorithm.
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
#' @param grid:        Matrix of the coordinates of grid points where to evaluate
#'                     the regression coefficients. dim(grid) = [n points] [2]
#' @param model:       Choose among ("midpoint","benchmark") or whichever other model to be tested
#' 
#' @return a six-element list with the following components:
#'         beta_c:     matrix of constant coefficients, evaluated on the grid points.
#'                     dim(beta_c) = [n constant coefficients] [n points]
#'         beta_s:     matrix of site-dependent coefficients, evaluated on the grid points.
#'                     dim(beta_s) = [n site-dependent coefficients] [n points]
#'         beta_e:     matrix of event-dependent coefficients, evaluated on the grid points.
#'                     dim(beta_e) = [n event-dependent coefficients] [n points]
#'

SEC_grid_creation = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, grid,
                             model){
  dist_e_sim = gw.dist(utm_ev_sp, grid, focus=0, p=2, theta=0, longlat=F)
  dist_s_sim = gw.dist(utm_st_sp, grid, focus=0, p=2, theta=0, longlat=F)
  
  N = length(y) #y vector of responses
  if (intercept == "c"){
    Xc = cbind(rep(1,N), Xc)
  } else if (intercept == "e"){
    Xe = cbind(rep(1,N), Xe)
  } else if (intercept == "s"){
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
  
  load(paste0(model,"/large_matrices/only_calibration_He.RData"))
  load(paste0(model,"/large_matrices/only_calibration_Hs.RData"))
  load(paste0(model,"/large_matrices/only_calibration_B.RData"))
  
  ## FUNCTIONS ----------------------------------------------------
  gauss_kernel = function(d, h){
    wgts = exp(-0.5*(d/h)^2)
    return(wgts)
  }
  
  my_fun_e = function(i, Xe, Hs, dist_e, bwe, gauss_kernel, y_tilde)
  {
    We = diag(gauss_kernel(dist_e[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    return(Ae %*% y_tilde)
  }
  
  my_fun_s = function(i, Xs, dist_s, bws, gauss_kernel, y_tilde_s)
  {
    Ws = diag(gauss_kernel(dist_s[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    return(As %*% y_tilde_s)
  }
  
  ## --------------------------------------------------------------
  
  #1)compute beta_c
  beta_c = solve(t(Xc)%*%t(B)%*%B%*%Xc)%*%t(Xc)%*%t(B)%*%B%*%y
  
  #2)compute y_tilde
  y_tilde = y-Xc%*%beta_c
  
  #3)compute beta_e for whole grid
  ncpu=4
  sfInit(par=TRUE,cp=ncpu)
  reps = 1:L
  (Start.Time <- Sys.time())
  beta_e = sfSapply(x = reps,
                fun = my_fun_e,
                Xe = Xe,
                Hs = Hs,
                dist_e = dist_e_sim,
                bwe = bwe,
                gauss_kernel = gauss_kernel,
                y_tilde = y_tilde)
  End.Time <- Sys.time()
  print(paste0("Building beta_e: ",round(End.Time - Start.Time, 2)))
  sfStop() #stop cluster parallelization
  
  #4)compute y_tilde_s
  y_tilde_s = (I-He)%*%y_tilde
  
  #5)compute beta_s for whole grid
  sfInit(par=TRUE,cp=ncpu)
  (Start.Time <- Sys.time())
  beta_s = sfSapply(x = reps,
                fun = my_fun_s,
                Xs = Xs,
                dist_s = dist_s_sim,
                bws = bws,
                gauss_kernel = gauss_kernel,
                y_tilde_s = y_tilde_s)
  End.Time <- Sys.time()
  print(paste0("Building beta_s: ",round(End.Time - Start.Time, 2)))
  sfStop() #stop cluster parallelization
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e)
  return(betas)
}