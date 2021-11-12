WMRFun = function(Gamma,se_Gamma,
                  alpha,se_alpha,
                  ldscore,R,MAF){
  
  alpha = alpha*sqrt(2*MAF(1-MAF))
  se_alpha=  se_alpha*sqrt(2*MAF(1-MAF))
  Gamma = Gamma*sqrt(2*MAF(1-MAF))
  se_Gamma = se_Gamma*sqrt(2*MAF(1-MAF))
  
  
  #initial estimate of beta
  
  beta_est = as.numeric(crossprod(Gamma,alpha)/crossprod(alpha))
  
  #step two: estimate tau
  diff_var = (se_Gamma^2+beta_est^2*se_alpha^2)
  chi_est = (Gamma-alpha*beta_est)^2/diff_var
  scale_ldscore = ldscore/diff_var
  tau_est = coefficients(lm(chi_est~scale_ldscore-1))
  tau_est = ifelse(tau_est>0,tau_est,0)
  
  
  
  #step three: estimate beta 
  W_inv = GetWtauMat(Gamma,se_Gamma,
                 alpha,se_alpha,
                 ldscore,tau_est,beta_est,R)
  alphaw_inv = t(solve(W_inv,alpha))
  awa = alphaw_inv%*%alpha
    
    #quadform(x= as.matrix(alpha),M = W)
  beta_est = awa^-1*
    crossprod(alphaw_inv,Gamma)
  #beta_se = sqrt(awa)
 # beta_var= (awa-quadform(x= as.matrix(se_alpha),M = W*R))^-1
  beta_var= (awa)^-1
  beta_se = sqrt(beta_var)
  #beta_se = sqrt(awa-quadform(x= as.matrix(se_alpha),M = diag(diag(W)))
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
  W_inv = diag(se_Gamma+beta_est*se_alpha)%*%R%*%diag(se_Gamma+beta_est*se_alpha)
  +tau*diag(ldscore)
  #W = solve(W_inv)
  #W[abs(W)<=1E-06] = 0
  return(W_inv)
}

