#'
#' Check of significance of the constant coefficients
#'
#' This function performs a non-parametric test to check whether a certain
#' regression coefficient, previously found to be constant, is equal to zero
#' or different from zero. In order to do so, two spatially varying regression
#' models are compared, one containing the coefficient to be checked-out and the
#' other one not containing it.
#' 
#' The formulation of the test is the following:
#' Let beta be a regression coefficient found to be non-spatially varying.
#' H0: beta  = 0
#' H1: beta != 0 
#'
#' @param coef_to_check:   name of the coefficient to be checked
#'                         (e.g. "intercept", "b1")
#' @param names_Xc         vector of names of the constant coefficients of the model              
#' @param Xc:              matrix of constant predictor variables
#' @param Xe:              matrix of event-varying predictor variables
#' @param Xs:              matrix of site-varying predictor variables
#' @param y:               response variable
#' @param bw1:             bandwidth for event or midpoint
#' @param bw2:             bandwidth for site
#' @param utm_1_sp:        utm coordinates of the events
#' @param utm_2_sp:        utm coordinates of the site
#' @param model:           choose among ("midpoint","benchmark") or whichever other model you're working with
#' 
#' @return p:              p-value of the non-parametric test
#' 



significance_check <- function(coef_to_check, names_Xc, Xc, Xe, Xs, y, bw1,
                               bw2, utm_1_sp, utm_2_sp, model)
{
  N = dim(Xc)[1]

  logic_vec = (names_Xc != coef_to_check)
  Xc = Xc[,logic_vec]
  
  if(coef_to_check=="intercept")
  {
    sec_null = SEC_calibration(Xc, Xe, Xs, y, "no", bw1, bw2, coordinates(utm_1_sp), coordinates(utm_2_sp), model, paste0("significance_",coef_to_check))
  } else {
    sec_null = SEC_calibration(Xc, Xe, Xs, y, "c", bw1, bw2, coordinates(utm_1_sp), coordinates(utm_2_sp), model, paste0("significance_",coef_to_check))
  }
  #compute R(H0)
  I= diag(1,N)
  B = sec_null$B
  if(coef_to_check=="intercept")
  {
    Xcc = Xc
  } else {
    Xcc = cbind(rep(1,N), Xc)
  }
  H0 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
  RH0 = t(I-H0)%*%(I-H0)
  epsilon= (I-H0)%*%y
  #load R(H1)
  load(paste0(model,"/large_matrices/RH1_stationary_rotD50pga.RData"))
  #compute T
  T0 = (t(y) %*% (RH0-RH1) %*% y) / (t(y) %*% RH1 %*% y)
  #permutations
  n_perm = 1000
  t_stat = rep(0,n_perm)
  set.seed(140996)
  print("Permutations")
  pb = progress_bar$new(total=n_perm, format = "  computing [:bar] :percent eta: :eta")
  for (i in 1:n_perm){
    eps_star = sample(epsilon)
    y_star = H0 %*% y + eps_star
    t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
    pb$tick()
  }
  p = sum(t_stat>as.numeric(T0))/n_perm
  
  return(p)
  
}