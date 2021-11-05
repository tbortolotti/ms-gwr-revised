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
#' @param emgwr:       list containing matrices He, Hs and B of the fully calibrated model
#' 
#' @return a six-element list with the following components:
#'         beta_c:     matrix of constant coefficients, evaluated on the grid points.
#'                     dim(beta_c) = [n constant coefficients] [n points]
#'         beta_s:     matrix of site-dependent coefficients, evaluated on the grid points.
#'                     dim(beta_s) = [n site-dependent coefficients] [n points]
#'         beta_e:     matrix of event-dependent coefficients, evaluated on the grid points.
#'                     dim(beta_e) = [n event-dependent coefficients] [n points]
#'         He:         matrix He
#'         Hs:         matrix Hs
#'         B:          matrix B
#'

SEC_grid_creation = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp, grid, emgwr){
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
  
  B = emgwr$B
  He = emgwr$He
  Hs = emgwr$Hs
  
  #1)compute beta_c
  beta_c = solve(t(Xc)%*%t(B)%*%B%*%Xc)%*%t(Xc)%*%t(B)%*%B%*%y
  
  #2)compute y_tilde
  y_tilde = y-Xc%*%beta_c
  
  #3)compute beta_e for whole grid
  print("Evaluation of event-dependent coefficients")
  pb = progress_bar$new(total=L, format = "  computing [:bar] :percent eta: :eta")
  for (i in 1:L){
    We = diag(gauss_kernel(dist_e_sim[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    beta_e[,i] = Ae %*% y_tilde
    #print(c("beta_e",i))
    pb$tick()
  }
  
  #4)compute y_tilde_s
  y_tilde_s = (I-He)%*%y_tilde
  
  #5)compute beta_s for whole grid
  print("Evaluation of site-dependent coefficients")
  pb = progress_bar$new(total=L, format = "  computing [:bar] :percent eta: :eta")
  for (i in 1:L){
    Ws = diag(gauss_kernel(dist_s_sim[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    beta_s[,i] = As %*% y_tilde_s
    #print(c("beta_s",i))
    pb$tick()
  }
  
  betas <- list("beta_c" = beta_c, 
                "beta_s" = beta_s,
                "beta_e" = beta_e,
                "He" = He,
                "Hs" = Hs,
                "B" = B)
  return(betas)
}