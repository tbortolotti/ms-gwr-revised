
#'
#' Evaluation of GCV statistic - constant intercept
#'
#' This function evaluates the gcv statistic for every site and event 
#' bandwidth combination.
#' Note that: the model has constant intercept
#'
#' @param bwe:         bandwidth for event
#' @param bws:         bandwidth for site
#' @param func:        one of the two functions: SEC_only_calibration/ESC_only_calibration
#' @param Xc:          matrix of constant predictor variables
#' @param Xe:          matrix of event-varying predictor variables
#' @param Xs:          matrix of site-varying predictor variables
#' @param y:           response variable
#' @param intercept:   "c" (constant)
#' @param utm_ev_sp:    utm coordinates of the events
#' @param utm_st_sp:   utm coordinates of the site   
#' 
#' @return a three-element list with the following components:
#'         bwe:    bandwidth for event, same as input
#'         bws:    bandwidth for site, same as input
#'         gcv:    value of the gcv statistic
#' 

gcv_mei_only_one = function(bwe, bws, func, Xc, Xe, Xs, y, intercept = "c", utm_ev_sp, utm_st_sp){
  temp = func(Xc, Xe, Xs, y, intercept, bwe, bws, utm_ev_sp, utm_st_sp)
  B = temp$B
  n = length(y)
  I = diag(1,n)
  Xcc = cbind(rep(1,n), Xc)
  H = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
  res = (I-H)%*%y
  gcv = sum(res^2/(rep(1, length(diag(H)))-diag(H))^2)
  
  solution = list("bwe" = bwe,
                  "bws" = bws,
                  "gcv" = gcv)
  return(solution)
}