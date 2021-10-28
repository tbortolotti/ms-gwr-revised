#'
#' Bandwidth selection via cross-validation
#'
#' This function performs the selection of the optimal bandwith of the gaussian kernel
#' used via cross-validation
#'
#' @param bw_min:      min value of the range of bandwidths
#' @param bw_max:      max value of the range of bandwidths
#' @param step:        step inside the range
#' @param f:           number of cross-validation folds
#' @param func:        one of the two functions: SEC_only_calibration/ESC_only_calibration
#' @param method:      either "sec" or "esc" (must be in accordance to func)
#' @param Xc, Xe, Xs:  predictor variables
#' @param y:           response variable
#' @param intercept:   either "c" (constant), "e" (event-dependent), "s" (site-dependent)
#' @param utm_ev_s:    utm coordinates of the events
#' @param utm_st_sp:   utm coordinates of the site   
#' 
#' @return a two-element list with the following components:
#'         best_bw:    value of the best selected bandwidth
#'         cvsum:      for each bandwidth, it is the sum, over all train-test partitions, of the RSS
#' 

bw_cv = function(bw_min, bw_max, step, f, func, method, Xc, Xe, Xs,
                 y,intercept, utm_ev_sp, utm_st_sp)
{
  bandwidths = seq(bw_min, bw_max, step)
  if (bandwidths[end(bandwidths)[1]]<bw_max-0.0001){
    bandwidths = c(bandwidths, bw_max)
  }
  cvss = matrix(0,length(bandwidths), f)
  batch = floor(length(y)/f)
  rest = length(y) - batch * f
  rest_const = rest
  for (w in 1:length(bandwidths)){
    i = bandwidths[w]
    print(c("Bandwidth", i))
    for (b in 1:f){
      if (rest >0){
        batch_t = batch + 1
        beginning = batch_t*(b-1)+1
        end = b*batch_t
        rest = rest - 1
      } else {
        beginning = rest_const + batch*(b-1) + 1 #(batch+1)*rest_const + batch*(b-1-rest_const)+1 
        end = rest_const + batch*b
      }
      #create subsets
      Xc_temp = as.matrix(Xc)[-(beginning:end),]
      Xe_temp = as.matrix(Xe)[-(beginning:end),]
      Xs_temp = as.matrix(Xs)[-(beginning:end),]
      y_temp = y[-(beginning:end)]
      coords_e_temp = utm_ev_sp[-(beginning:end),]
      coords_s_temp = utm_st_sp[-(beginning:end),]
      
      Xc_test = as.matrix(Xc)[(beginning:end),]
      Xe_test = as.matrix(Xe)[(beginning:end),]
      Xs_test = as.matrix(Xs)[(beginning:end),]
      y_test = y[(beginning:end)]
      coords_e_test = utm_ev_sp[(beginning:end),]
      coords_s_test = utm_st_sp[(beginning:end),]
      
      if (batch == 1){
        pcoords = c(coords_e_test, coords_s_test)
      } else {
        pcoords = cbind(coords_e_test, coords_s_test)
      }
      
      temp = func(Xc_temp, Xe_temp, Xs_temp, y_temp, intercept,
                  i, i, coords_e_temp, coords_s_temp)
      if (method == "sec"){
        prediction = emgwr_prediction_points_no_param(Xc_temp, Xe_temp, Xs_temp, y_temp, temp, "sec", i, i,
                                                      coords_e_temp, coords_s_temp, Xc_test, Xe_test, Xs_test,
                                                      pcoords, 0.05)
      } else {
        prediction = emgwr_prediction_points_no_param(Xc_temp, Xe_temp, Xs_temp, y_temp, temp, "esc", i, i,
                                                      coords_e_temp, coords_s_temp, Xc_test, Xe_test, Xs_test,
                                                      pcoords, 0.05)
      }
      
      
      
      #compute residuals
      y_res = rep(0,length(y_test))
      y_hat = prediction$y0
      y_res = (y_hat - y_test)
      cvss[w,b] = Norm(y_res,2)^2
    }
  }
  
  cvsum = matrix(0, length(bandwidths), 2)
  for (w in 1:length(bandwidths)){
    i = bandwidths[w]
    cvsum[w,] = c(i,sum(cvss[w,]))
  }
  best = min(cvsum[,2])
  best_bw = cvsum[(cvsum[,2]==best),1]
  
  solution = list("best_bw" = best_bw, 
                  "cvsum" = cvsum
  )
  return(solution)
}
