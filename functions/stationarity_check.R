stationarity_check <- function(coef_to_check, regs, y, bwe, bws, utm_ev_sp, utm_st_sp)
{
  ## coefs
  b1 = regs$b1
  b2 = regs$b2
  c1 = regs$c1
  c2 = regs$c2
  c3 = regs$c3
  f1 = regs$f1
  f2 = regs$f2
  k  = regs$k
  
  e_dependent = cbind(c2,c3)
  s_dependent = cbind(b1,b2,f1,f2,c1,k)
  
  names_e_dependent = c("c2","c3")
  logic_vec = (names_e_dependent != coef_to_check)
  e_dependent = e_dependent[,logic_vec]
  
  names_s_dependent = c("b1","b2","f1","f2","c1","k")
  logic_vec = (names_s_dependent != coef_to_check)
  s_dependent = s_dependent[,logic_vec]
  
  Xc = coef_to_check
  Xe = e_dependent
  Xs = s_dependent
  
  sec = SEC_only_calibration(Xc, Xe, Xs, y, "c", bwe, bws, coordinates(utm_ev_sp), coordinates(utm_st_sp))
  #compute R(H0)
  I= diag(1,n_sample)
  B = sec$B
  Xcc = cbind(rep(1,n_sample), Xc)
  H0 = I - B + B %*% Xcc %*% solve(t(Xcc)%*%t(B)%*%B%*%Xcc) %*% t(Xcc) %*% t(B)%*% B
  RH0 = t(I-H0)%*%(I-H0)
  epsilon= (I-H0)%*%y
  #load R(H1)
  RH1 = readRDS("RH1_only_intercept_rotD50pga.RData")
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