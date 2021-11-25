WMRFun = function(Gamma,se_Gamma,
                  alpha,se_alpha,
                  ldscore,R,MAF){
  #initial estimate of beta
  alpha = alpha*sqrt(2*MAF*(1-MAF))
  se_alpha=  se_alpha*sqrt(2*MAF*(1-MAF))
  Gamma = Gamma*sqrt(2*MAF*(1-MAF))
  se_Gamma = se_Gamma*sqrt(2*MAF*(1-MAF))
  beta_est = as.numeric(crossprod(Gamma,alpha)/crossprod(alpha))
  
  #step two: estimate tau
  diff_var = (se_Gamma^2+beta_est^2*se_alpha^2)
  chi_est = (Gamma-alpha*beta_est)^2/diff_var
  scale_ldscore = ldscore/diff_var
  tau_est = coefficients(lm(chi_est~scale_ldscore-1))
  tau_est = ifelse(tau_est>0,tau_est,0)

  beta_var= (awa-quadform(x= as.matrix(se_alpha),M = W*R))^-1
  #beta_var= awa^-1
  beta_se = sqrt(beta_var)
  # best_est = (alpha%*%W%*%Gamma)/(alpha%*%W%*%alpha)
  # print(best_est)
  # print(beta_est)
  return(c(beta_est,beta_se,beta_low=beta_est-1.96*beta_se,
           beta_high = beta_est+1.96*beta_se))
}

GetWtauMat = function(Gamma,se_Gamma,
                      alpha,se_alpha,
                      ldscore,tau,beta_est,R){
  #  W = diag(1/(se_Gamma^2+beta_est^2*se_alpha^2+tau*ldscore))
  W_inv = diag(se_Gamma)%*%R%*%diag(se_Gamma)+
    beta_est^2*diag(se_alpha)%*%R%*%diag(se_alpha)
  
  #W_inv = (se_Gamma+beta_est*se_alpha)%*%(se_Gamma+beta_est*se_alpha)
  #+tau*diag(ldscore)
  W = solve(W_inv)
  #W[abs(W)<=1E-06] = 0
  return(W)
}


ObjFun <- function(Gamma,alpha,W,beta_est){
  return(quadform(x= as.matrix(Gamma-beta_est*alpha),M = W))
}