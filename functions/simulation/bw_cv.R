#'
#' Bandwidth selection via cross-validation
#'
#' This function performs a cross-validation to choose the optimal bandwidths, which are those providing the
#' minimum of the mean residual sum of squares across all partitions.
#' 
#' @param bw_min:      minimum value of the grid of bandwidths to try
#' @param bw_max:      maximum value of the grid of bandwidths to try
#' @param step:        step of the bandwidths grid
#' @param f            number of batches of the f-fold cross-validation
#' @param func:        either ESC_calibration or SEC_calibration, depending on which calibration algorithm is
#'                     to be used
#' @param method:      either "esc" or "sec", in accordance to what is specified for func
#' @param Xc:          matrix of constant predictor variables
#' @param Xe:          matrix of event-varying predictor variables
#' @param Xs:          matrix of site-varying predictor variables
#' @param y:           response variable
#' @param intercept:   either "c" (constant), "e" (event-dependent), "s" (site-dependent)
#' @param utm_1_sp:    utm coordinates of the events or midpoints
#' @param utm_2_sp:    utm coordinates of the site
#'                      
#' @return a two-element list with the following components:
#'         best_bw:    the optimal bandwidths for site and event
#'         cvsum:      value of the mean residual sum of squares, resulting from CV
#'


bw_cv = function(bw_min, bw_max, step, f, func, method, Xc, Xe, Xs, y,intercept, utm_1_sp, utm_2_sp){
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
    pb = progress_bar$new(total=f, format = "  computing [:bar] :percent eta: :eta")
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
      coords_e_temp = utm_1_sp[-(beginning:end),]
      coords_s_temp = utm_2_sp[-(beginning:end),]
      
      Xc_test = as.matrix(Xc)[(beginning:end),]
      Xe_test = as.matrix(Xe)[(beginning:end),]
      Xs_test = as.matrix(Xs)[(beginning:end),]
      y_test = y[(beginning:end)]
      coords_e_test = utm_1_sp[(beginning:end),]
      coords_s_test = utm_2_sp[(beginning:end),]
      
      if (batch == 1){
        pcoords = c(coords_e_test, coords_s_test)
      } else {
        pcoords = cbind(coords_e_test, coords_s_test)
      }
      
      temp = func(Xc_temp, Xe_temp, Xs_temp, y_temp, intercept, i, i, coords_e_temp, coords_s_temp)
      if (method == "sec"){
        prediction = emgwr_calibration_prediction(Xc_temp, Xe_temp, Xs_temp, y_temp, temp, "sec", i, i,
                                                      coords_e_temp, coords_s_temp, Xc_test, Xe_test, Xs_test,
                                                      pcoords, 0.05)
      } else {
        prediction = emgwr_calibration_prediction(Xc_temp, Xe_temp, Xs_temp, y_temp, temp, "esc", i, i,
                                                      coords_e_temp, coords_s_temp, Xc_test, Xe_test, Xs_test,
                                                      pcoords, 0.05)
      }
      
      
      
      #compute residuals
      y_res = rep(0,length(y_test))
      y_hat = prediction$y0
      y_res = (y_hat - y_test)
      cvss[w,b] = Norm(y_res,2)^2
      
      pb$tick()
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
