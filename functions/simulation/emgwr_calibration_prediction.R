#'
#' Model calibration and prediction via MS-GWR, either via SEC or ESC algorithms
#'
#' This function builds matrices He, Hs and B that are used for the MS-GWR
#' via SEC algorithm or ESC algorithm, and evaluates the predictions for points the covariates
#' and coordinates of which are given as input
#'
#' @param Xc:          matrix of constant predictor variables on which the model is calibrated (train)
#' @param Xe:          matrix of event-varying predictor variables on which the model is calibrated (train)
#' @param Xs:          matrix of site-varying predictor variables on which the model is calibrated (train)
#' @param y:           response variable on which the model is calibrated (train)
#' @param emgwr:       model resulting from the performance of either SEC_calibration or ESC_calibration, according to the
#'                     algorithm of calibration and prediction
#' @param method:      either "sec" or "esc" according to what is specified for emgwr
#' @param bw1:         bandwidth for event or midpoint
#' @param bw2:         bandwidth for site
#' @param utm_1_sp:    utm coordinates of the events or midpoints of points on which the model is calibrated (train)
#' @param utm_2_sp:    utm coordinates of the site of points on which the model is calibrated (train)
#' @param pc:          matrix (or vector) of constant predictor variables for which the model predicts (test)
#' @param pe:          matrix (or vector) of event-dependent predictor variables for which the model predicts (test)
#' @param ps:          matrix (or vector) of site-dependent predictor variables for which the model predicts (test)
#' @param pcoords:     coordinates of the points for which we want the prediction
#' @param alpha:        level of the prediction interval associated to the point-wise predictions
#' 
#' @return a six-element list with the following components:
#'         s0:            s0
#'         y0:            vector of fitted values
#'         betas:         vector of estimated regression coefficients
#'         interval:      prediction interval of level 1-alpha
#'         var0:          estimate of Var(y0), i.e. the prediction uncertainty
#'         sigma2hat:     estimate of the variability of the residuals
#'


emgwr_calibration_prediction = function(Xc, Xe, Xs, y, emgwr, method, bw1, bw2, utm_1_sp, utm_2_sp, pc, pe, ps, pcoords, alpha){
  
  N = length(y)
  K = dim(pcoords)[1]
  if (is.null(K)){
    K = 1
  }
  Xc = cbind(rep(1,N), Xc)
  if (K==1 | ((is.null(dim(pc)) & length(pc)!= K & length(pe)!= K & length(ps)!= K))){ #only one point to predict
    x0 = matrix(c(1, pc, pe, ps), ncol = 1)
  } else {
    x0 = cbind(1, pc, pe, ps)
  }
  
  I = diag(rep(1,N))
  
  Xc = as.matrix(Xc)
  n_c = dim(Xc)[2] #number constant of covariates
  Xe = as.matrix(Xe)
  n_e = dim(Xe)[2] #number of event-dependent covariates
  Xs = as.matrix(Xs)
  n_s = dim(Xs)[2] #number of site-dependent covariates
  
  He = emgwr$He
  Hs = emgwr$Hs
  B = emgwr$B
  
  H_hat = I - B + B %*% Xc %*% solve(t(Xc)%*%t(B)%*%B%*%Xc) %*% t(Xc) %*% t(B)%*% B
  res = (I-H_hat)%*%y
  delta1 = N-2*tr(H_hat)+tr(t(H_hat)%*%H_hat)
  delta2 = tr((t(I-H_hat)%*%(I-H_hat)) %*% (t(I-H_hat)%*%(I-H_hat)))
  rss = sum(res^2)
  sigma2hat = rss/delta1
  sigmahat = sqrt(sigma2hat)
  
  #1)compute beta_c
  Ac = solve(t(Xc)%*%t(B)%*%B%*%Xc)%*%t(Xc)%*%t(B)%*%B
  beta_c = Ac%*%y
  
  #2)compute y_tilde
  y_tilde = y-Xc%*%beta_c
  
  if (K==1){
    dist_e_sim = gw.dist(utm_1_sp, t(as.matrix(pcoords)[1:2]), focus=0, p=2, theta=0, longlat=F)
    dist_s_sim = gw.dist(utm_2_sp, t(as.matrix(pcoords)[3:4]), focus=0, p=2, theta=0, longlat=F)
    
    if (method == "esc") {
      Ws = diag(gauss_kernel(dist_s_sim[,1],bw2))
      As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
      beta_s = As %*% y_tilde
      
      y_tilde_s = (I-Hs)%*%y_tilde
      
      We = diag(gauss_kernel(dist_e_sim[,1],bw1))
      Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
      beta_e = Ae %*% y_tilde_s
      
      Q = rbind(Ac, Ae%*%(I-Hs)%*%(I-Xc%*%Ac), As%*%(I-Xc%*%Ac))
      
      s0 = t(x0)%*%Q%*%t(Q)%*%x0
      
      betas = c(beta_c, beta_e, beta_s)
      y0 = t(x0)%*%betas
      t_alpha2 = qt(1-alpha/2, delta1^2/delta2)
      var0 = sigma2hat*s0
      interval = c(y0-sigmahat*sqrt(s0+1)*t_alpha2, y0, y0+sigmahat*sqrt(s0+1)*t_alpha2)
      
    } else {
      We = diag(gauss_kernel(dist_e_sim[,1],bw1))
      Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
      beta_e = Ae %*% y_tilde
      
      y_tilde_s = (I-He)%*%y_tilde
      
      Ws = diag(gauss_kernel(dist_s_sim[,1],bw2))
      As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
      beta_s = As %*% y_tilde_s
      
      Q = rbind(Ac, Ae%*%(I-Xc%*%Ac), As%*%(I-He)%*%(I-Xc%*%Ac))
      
      s0 = t(x0)%*%Q%*%t(Q)%*%x0
      
      betas = c(beta_c, beta_e, beta_s)
      y0 = t(x0)%*%betas
      t_alpha2 = qt(1-alpha/2, delta1^2/delta2)
      var0 = sigma2hat*s0
      interval = c(y0-sigmahat*sqrt(s0+1)*t_alpha2, y0, y0+sigmahat*sqrt(s0+1)*t_alpha2)
    }
  }
  
  else if (K > 1){
    betas = matrix(0,K,n_c+n_e+n_s)
    y0 = rep(0,K)
    var0 = rep(0,K)
    s0 = rep(0,K)
    interval = matrix(0,K,3)
    
    if (method == "esc"){
      for (k in 1:K){
        dist_e_sim = gw.dist(utm_1_sp, pcoords[k,1:2], focus=0, p=2, theta=0, longlat=F)
        dist_s_sim = gw.dist(utm_2_sp, pcoords[k,3:4], focus=0, p=2, theta=0, longlat=F)
        
        Ws = diag(gauss_kernel(dist_s_sim[,1],bw2))
        As = (solve((t(Xs)%*%(t((I-He)))%*%Ws%*%(I-He)%*%Xs))) %*% t(Xs) %*% (t((I-He))) %*% Ws %*% (I-He)
        beta_s = As %*% y_tilde
        
        y_tilde_s = (I-Hs)%*%y_tilde
        
        We = diag(gauss_kernel(dist_e_sim[,1],bw1))
        Ae = (solve((t(Xe)%*%We%*%Xe))) %*% t(Xe) %*% We
        beta_e = Ae %*% y_tilde_s
        
        Q = rbind(Ac, Ae%*%(I-Hs)%*%(I-Xc%*%Ac),As%*%(I-Xc%*%Ac))
        
        s0[k] = t(x0[k,])%*%Q%*%t(Q)%*%x0[k,]
        
        betas[k,] = c(beta_c, beta_e, beta_s)
        y0[k] = x0[k,]%*%betas[k,]
        t_alpha2 = qt(1-alpha/2, delta1^2/delta2)
        var0[k] = sigma2hat*s0[k]
        interval[k,] = c(y0[k]-sigmahat*sqrt(s0[k]+1)*t_alpha2, y0[k], y0[k]+sigmahat*sqrt(s0[k]+1)*t_alpha2)
      }
    } else {
      for (k in 1:K){
        dist_e_sim = gw.dist(utm_1_sp, pcoords[k,1:2], focus=0, p=2, theta=0, longlat=F)
        dist_s_sim = gw.dist(utm_2_sp, pcoords[k,3:4], focus=0, p=2, theta=0, longlat=F)
        
        We = diag(gauss_kernel(dist_e_sim[,1],bw1))
        Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
        beta_e = Ae %*% y_tilde
        
        y_tilde_s = (I-He)%*%y_tilde
        
        Ws = diag(gauss_kernel(dist_s_sim[,1],bw2))
        As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
        beta_s = As %*% y_tilde_s
        
        Q = rbind(Ac, Ae%*%(I-Xc%*%Ac), As%*%(I-He)%*%(I-Xc%*%Ac))
        s0[k] = t(x0[k,])%*%Q%*%t(Q)%*%x0[k,]
        
        betas[k,] = c(beta_c, beta_e, beta_s)
        y0[k] = x0[k,]%*%betas[k,]
        t_alpha2 = qt(1-alpha/2, delta1^2/delta2)
        var0[k] = sigma2hat*s0[k]
        interval[k,] = c(y0[k]-sigmahat*sqrt(s0[k]+1)*t_alpha2, y0[k], y0[k]+sigmahat*sqrt(s0[k]+1)*t_alpha2)
      }
    }
  }
  
  
  result <- list("s0" = s0, 
                 "y0" = y0,
                 "betas" = betas,
                 "interval" = interval,
                 "var0" = var0,
                 "sigma2hat" = sigma2hat)
  return(result)
}
