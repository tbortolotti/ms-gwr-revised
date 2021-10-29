significance_check <- function(coef_to_check, names_Xc, Xc, Xe, Xs, y, bwe,
                               bws, utm_ev_sp, utm_st_sp)
{

  names_Xc = c("b1","b2","f1","f2","c1")
  logic_vec = (names_Xc != coef_to_check)
  Xc = Xc[,logic_vec]
  
  if(coef_to_check=="intercept")
  {
    sec_null = SEC_no_intercept_calibration(Xc, Xe, Xs, y, bwe, bws, coordinates(utm_ev_sp), coordinates(utm_st_sp))
  } else {
    sec_null = SEC_only_calibration(Xc, Xe, Xs, y, "c", bwe, bws, coordinates(utm_ev_sp), coordinates(utm_st_sp))
  }
  #compute R(H0)
  I= diag(1,n_sample)
  B = sec_null$B
  if(coef_to_check=="intercept")
  {
    Xcc = Xc
  } else {
    Xcc = cbind(rep(1,n_sample), Xc)
  }
  H0 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
  RH0 = t(I-H0)%*%(I-H0)
  epsilon= (I-H0)%*%y
  #load R(H1)
  RH1 = readRDS("RH1_stationary_rotD50pga.RData")
  #compute T
  T0 = (t(y) %*% (RH0-RH1) %*% y) / (t(y) %*% RH1 %*% y)
  #permutations
  n_perm = 1000
  t_stat = rep(0,n_perm)
  print("Permutations")
  pb = progress_bar$new(total=n_perm, format = "  computing [:bar] :percent eta: :eta")
  pb$tick(0)
  for (i in 1:n_perm){
    eps_star = sample(epsilon)
    y_star = H0 %*% y + eps_star
    t_stat[i] = (t(y_star) %*% (RH0-RH1) %*% y_star) / (t(y_star) %*% RH1 %*% y_star)
    pb$tick()
  }
  p = sum(t_stat>as.numeric(T0))/n_perm
  
  return(p)
  
}