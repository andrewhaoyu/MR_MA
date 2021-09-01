WMRFun = function(Gamma,se_Gamma,
                  alpha,se_alpha,
                  ldscore,R){
  #initial estimate of beta
  beta_est = crossprod(Gamma,alpha)/crossprod(alpha)
  
  #step two: estimate tau
  diff_var = (se_Gamma^2+beta_est^2*se_alpha^2)
  chi_est = (Gamma-alpha*beta_est)^2/diff_var
  scale_ldscore = ldscore/diff_var
  tau_est = coefficients(lm(chi_est~scale_ldscore-1))
  
}
# GetWtauMat = function(Gamma,se_Gamma,
#                       alpha,se_alpha,
#                       ldscore,tau,beta_est){
#   W = 1/diag(se_Gamma^2+beta_est^2*se_alpha^2+tau*ldscore)
#   return(W)
# }