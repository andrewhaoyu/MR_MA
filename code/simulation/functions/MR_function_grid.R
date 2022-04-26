WMRFun = function(Gamma,se_Gamma,
                  alpha,se_alpha,
                  ldscore,R,MAF,tau){
  #initial estimate of beta
  alpha = alpha*sqrt(2*MAF*(1-MAF))
  se_alpha=  se_alpha*sqrt(2*MAF*(1-MAF))
  Gamma = Gamma*sqrt(2*MAF*(1-MAF))
  se_Gamma = se_Gamma*sqrt(2*MAF*(1-MAF))
  beta_est = as.numeric(crossprod(Gamma,alpha)/crossprod(alpha))
  
  beta_old = beta_est  
  #step two: estimate tau
  diff_var = (se_Gamma^2+beta_est^2*se_alpha^2)
  chi_est = (Gamma-alpha*beta_est)^2/diff_var
  scale_ldscore = ldscore/diff_var
  tau_est = coefficients(lm(chi_est~scale_ldscore))[2]
  tau_est = ifelse(tau_est>1E-05,tau_est,1E-05)
  tau = tau_est
  #tau_est = tau
  library(cPCG)
  
  
  
  SGRSG = GetSGRSG(se_Gamma,R)
  SARSA = GetSARSA(se_Gamma,R)
  TauL = GetTauL(tau_est,ldscore)
  
  
  
  V = GetVmat(SGRSG,SARSA,
              TauL,beta_est)
  #alpha %*%W as wa 
  wa = cgsolve(A=V, alpha)
  #t(alpha) %*%W %*% alpha as awa
  awa = crossprod(alpha,wa)
  #t(alpha) %*%W %*% gamma as awg
  awg = crossprod(Gamma,wa)
  beta_est = as.numeric(awg/awa)
  
  beta_new = beta_est
  error = abs(beta_new-beta_old)/abs(beta_old)
  print(error)
  while(error>0.05){
    beta_old = beta_est
    #step two: estimate tau
    diff_var = (se_Gamma^2+beta_est^2*se_alpha^2)
    chi_est = (Gamma-alpha*beta_est)^2/diff_var
    scale_ldscore = ldscore/diff_var
    tau_est = coefficients(lm(chi_est~scale_ldscore))[2]
    tau_est = ifelse(tau_est>1E-05,tau_est,1E-05)
    TauL = GetTauL(tau_est,ldscore)
    V = GetVmat(SGRSG,SARSA,
                TauL,beta_est)
    wa = cgsolve(A=V, alpha)
    #t(alpha) %*%W %*% alpha as awa
    awa = crossprod(alpha,wa)
    #t(alpha) %*%W %*% gamma as awg
    awg = crossprod(Gamma,wa)
    beta_est = as.numeric(awg/awa)
    beta_new = beta_est
    error = abs(beta_new-beta_old)/abs(beta_old)
    print(error)
  }
  
  #se_alpha %*% W as wse
  wse = cgsolve(A=V, b = se_alpha)
  beta_var= (awa-crossprod(se_alpha,wse))^-1
  beta_se = sqrt(beta_var)
  # best_est = (alpha%*%W%*%Gamma)/(alpha%*%W%*%alpha)
  # print(best_est)
  # print(beta_est)
  return(c(beta_est,beta_se,beta_low=beta_est-1.96*beta_se,
           beta_high = beta_est+1.96*beta_se))
}








GetSGRSG = function(se_Gamma,R){
  return(t(se_Gamma*t(se_Gamma*R)))
}
GetSARSA = function(se_Gamma,R){
  return(t(se_alpha*t(se_alpha*R)))
}
GetTauL = function(tau,ldscore){
  return(tau*diag(ldscore))
}

GetVmat = function(SGRSG,SARSA,
                   TauL,beta_est){
  #  W = diag(1/(se_Gamma^2+beta_est^2*se_alpha^2+tau*ldscore))
  # W_inv = diag(se_Gamma)%*%R%*%diag(se_Gamma)+
  #   beta_est^2*diag(se_alpha)%*%R%*%diag(se_alpha)
 V = SGRSG+
    as.numeric(beta_est)^2*SARSA+TauL
  
  
  
  # W_inv =  diag(se_Gamma)%*%R%*%diag(se_Gamma)+
  #   beta_est^2*diag(se_alpha)%*%R%*%diag(se_alpha)+tau*diag(ldscore)
  #W_inv = (se_Gamma+beta_est*se_alpha)%*%(se_Gamma+beta_est*se_alpha)
  #+tau*diag(ldscore)
  #W[abs(W)<=1E-06] = 0
  return(V)
}


ObjFun <- function(Gamma,alpha,W,beta_est){
  return(quadform(x= as.matrix(Gamma-beta_est*alpha),M = W))
}










# WMRFun = function(Gamma,se_Gamma,
#                   alpha,se_alpha,
#                   ldscore,R,MAF,tau){
#   #initial estimate of beta
#   alpha = alpha*sqrt(2*MAF*(1-MAF))
#   se_alpha=  se_alpha*sqrt(2*MAF*(1-MAF))
#   Gamma = Gamma*sqrt(2*MAF*(1-MAF))
#   se_Gamma = se_Gamma*sqrt(2*MAF*(1-MAF))
#   beta_est = as.numeric(crossprod(Gamma,alpha)/crossprod(alpha))
#   beta_old = beta_est  
#   #step two: estimate tau
#   diff_var = (se_Gamma^2+beta_est^2*se_alpha^2)
#   chi_est = (Gamma-alpha*beta_est)^2/diff_var
#   scale_ldscore = ldscore/diff_var
#   tau_est = coefficients(lm(chi_est~scale_ldscore-1))
#   tau_est = ifelse(tau_est>0,tau_est,0)
#   
#   
#   #tau_est = tau
#   W = GetWtauMat(Gamma,se_Gamma,
#                  alpha,se_alpha,
#                  ldscore,tau_est,beta_est,R)
#   awa = quadform(x= as.matrix(alpha),M = W)
#   beta_est = as.numeric(awa^-1*
#                           crossprod(t(crossprod(alpha,W)),Gamma))
#   beta_new = beta_est
#   error = abs(beta_new-beta_old)/abs(beta_old)
#   print(error)
#   while(error>0.05){
#     beta_old = beta_est
#     diff_var = (se_Gamma^2+beta_est^2*se_alpha^2)
#     chi_est = (Gamma-alpha*beta_est)^2/diff_var
#     scale_ldscore = ldscore/diff_var
#     tau_est = coefficients(lm(chi_est~scale_ldscore-1))
#     tau_est = ifelse(tau_est>0,tau_est,0)
#     #tau_est = tau
#     W = GetWtauMat(Gamma,se_Gamma,
#                    alpha,se_alpha,
#                    ldscore,tau_est,beta_est,R)
#     awa = quadform(x= as.matrix(alpha),M = W)
#     beta_est = as.numeric(awa^-1*
#                             crossprod(t(crossprod(alpha,W)),Gamma))
#     beta_new = beta_est
#     error = abs(beta_new-beta_old)/abs(beta_old)
#     print(error)
#   }
#   #beta_var= (awa-quadform(x= as.matrix(se_alpha),M = W*R))^-1
#   beta_var= (awa-quadform(x= as.matrix(se_alpha),M = W))^-1
#   #beta_var= awa^-1
#   beta_se = sqrt(beta_var)
#   
#   return(c(beta_est,beta_se,beta_low=beta_est-1.96*beta_se,
#            beta_high = beta_est+1.96*beta_se))
# }


# GetWtauMat = function(Gamma,se_Gamma,
#                       alpha,se_alpha,
#                       ldscore,tau,beta_est,R){
#   #  W = diag(1/(se_Gamma^2+beta_est^2*se_alpha^2+tau*ldscore))
#   # W_inv = diag(se_Gamma)%*%R%*%diag(se_Gamma)+
#   #   beta_est^2*diag(se_alpha)%*%R%*%diag(se_alpha)
#   # W_inv = t(se_Gamma*t(se_Gamma*R))+
#   #   beta_est^2*t(se_alpha*t(se_alpha*R))+tau*diag(ldscore)
#   #   
#   
#   
#   W_inv =  diag(se_Gamma)%*%R%*%diag(se_Gamma)+
#     beta_est^2*diag(se_alpha)%*%R%*%diag(se_alpha)+tau*diag(ldscore)
#   #W_inv = (se_Gamma+beta_est*se_alpha)%*%(se_Gamma+beta_est*se_alpha)
#   #+tau*diag(ldscore)
#   W = solve(W_inv)
#   #W[abs(W)<=1E-06] = 0
#   return(W)
# }
# 
