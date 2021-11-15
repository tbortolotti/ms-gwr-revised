
#'
#' Evaluation of GCV statistic - constant intercept
#'
#' This function evaluates the gcv statistic for every site and event 
#' bandwidth combination.
#' Note that: the model has constant intercept
#'
#' @param bw1:         bandwidth for event or midpoint
#' @param bw2:         bandwidth for site
#' @param func:        one of the two functions: SEC_only_calibration/ESC_only_calibration
#' @param Xc:          matrix of constant predictor variables
#' @param Xe:          matrix of event-varying predictor variables
#' @param Xs:          matrix of site-varying predictor variables
#' @param y:           response variable
#' @param intercept:   "c" (constant)
#' @param utm_1_sp:    utm coordinates of the events or midpoints
#' @param utm_2_sp:    utm coordinates of the site   
#' 
#' @return a three-element list with the following components:
#'         bw1:    bandwidth for event, same as input
#'         bw2:    bandwidth for site, same as input
#'         gcv:    value of the gcv statistic
#' 

gcv_mei_only_one = function(bw1, bw2, func, Xc, Xe, Xs, y, intercept = "c", utm_1_sp, utm_2_sp, model){
  temp = func(Xc, Xe, Xs, y, intercept, bw1, bw2, utm_1_sp, utm_2_sp, model, "bw_selection")
  B = temp$B
  n = length(y)
  I = diag(1,n)
  Xcc = cbind(rep(1,n), Xc)
  H = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
  res = (I-H)%*%y
  gcv = sum(res^2/(rep(1, length(diag(H)))-diag(H))^2)
  
  solution = list("bw1" = bw1,
                  "bw2" = bw2,
                  "gcv" = gcv)
  return(solution)
}