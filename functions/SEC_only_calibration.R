#SEC_only_calibration----
SEC_only_calibration = function(Xc, Xe, Xs, y,intercept, bwe, bws, utm_ev_sp, utm_st_sp){
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
  
  #create Hs
  print("Building Hs")
  pb = progress_bar$new(total=N, format = "  computing [:bar] :percent eta: :eta")
  for (i in 1:N){
    Ws = diag(gauss_kernel(dist_s_sim_cal[,i],bws))
    As = (solve((t(Xs)%*%Ws%*%Xs))) %*% t(Xs) %*% Ws
    Hs[i,] = Xs[i,] %*% As
    pb$tick()
  }
  
  #create He
  print("Building He")
  pb = progress_bar$new(total=N, format = "  computing [:bar] :percent eta: :eta")
  for (i in 1:N){
    We = diag(gauss_kernel(dist_e_sim_cal[,i],bwe))
    Ae = (solve((t(Xe)%*%(t((I-Hs)))%*%We%*%(I-Hs)%*%Xe))) %*% t(Xe) %*% (t((I-Hs))) %*% We %*% (I-Hs)
    He[i,] = Xe[i,] %*% Ae
    pb$tick()
  }
  
  #create B
  B = I - He - Hs + Hs %*% He
  
  calibration <- list("He" = He,
                      "Hs" = Hs,
                      "B" = B)
  
  return(calibration)
}