eff_alpha = 0.02
beta = 0.15
N = 60000
sigma_alpha = 1/N
sigma_Gamma = 1/N

n.rep = 1000

n.iv = 30

sigma_tau = 0.01^2

beta_est = rep(0,n.rep)
beta_cover = rep(0,n.rep)
beta_se = rep(0,n.rep)


beta_est_Raps = rep(0,n.rep)
beta_cover_Raps = rep(0,n.rep)

beta_est_IVW = rep(0,n.rep)
beta_cover_IVW = rep(0,n.rep)

beta_est_raps = rep(0,n.rep)
beta_cover_raps = rep(0,n.rep)


library(mr.raps)
library(MASS)
library(Rfast)

#independent IVs, no pleotropic
for(k in 1:n.rep){
  if(k%%100){print(k)}
  R = diag(n.iv)
  #R[upper.tri(R)] = R[lower.tri(R)] = 0.5
  #mean_alpha = c(eff_alpha,rep(eff_alpha*0.5,n.iv-1))
  mean_alpha = rep(eff_alpha,n.iv)
  
  
  
  alpha = as.vector(rmvnorm(1,mu = mean_alpha,sigma = R*sigma_alpha))
  #alpha[1] = rnorm(1,mean_alpha,sd = sqrt(sigma_alpha))
  var_alpha = rep(sigma_alpha,n.iv)
  alpha_temp = rmvnorm(1,mu = mean_alpha,sigma = R*sigma_alpha)
  #alpha_temp[1] = rnorm(1,mean_alpha,sd = sqrt(sigma_alpha))
  error_y = rmvnorm(1,mu = rep(0,n.iv),sigma = R*sigma_Gamma)
  Gamma = as.vector(alpha_temp*beta  + error_y)
  var_Gamma = rep(sigma_Gamma,n.iv)
  #R = diag(n.iv)
  MR_result <- MRWeight(Gamma,
                        var_Gamma,
                        alpha,
                        var_alpha,R)
  beta_est[k] = MR_result[[1]]
  beta_cover[k] = ifelse(MR_result[[2]]<=beta&MR_result[[3]]>=beta,1,0)
  beta_se[k] = MR_result[[4]]
  
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                      alpha,var_alpha)
  beta_est_IVW[k] = IVW_c_temp[1]
  beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
                            IVW_c_temp[3]>=beta,1,0)
  
  # raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
  #                                          beta.outcome = Gamma,
  #                                          se.exposure = sqrt(var_alpha),
  #                                          se.outcome = sqrt(var_Gamma)),
  #                        diagnostics = F)
  # 
  # beta_est_raps[k] = raps_est = raps_result$beta.hat
  # raps_se = raps_result$beta.se
  # beta_cover_raps[k] = ifelse(raps_est-1.96*raps_se<=beta&
  #                               raps_est+1.96*raps_se>=beta,1,0)
 
  # 
  # 
  raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
                                           beta.outcome = Gamma,
                                           se.exposure = sqrt(var_alpha),
                                           se.outcome = sqrt(var_Gamma)),
                         diagnostics = F)

  beta_est_Raps[k] = raps_result$beta.hat
  beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                             raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
 
}
# cov(beta_est,beta_cover)
# 
# 0.5*(sigma_Gamma+0.15^2*sigma_alpha)


paste0(round(mean(beta_est),3)," (",round(sd(beta_est),3),")")
mean(beta_cover)
mean(beta_se)
sd(beta_est)

paste0(round(mean(beta_est_IVW),3)," (",round(sd(beta_est_IVW),3),")")
mean(beta_cover_IVW)
paste0(round(mean(beta_est_Raps),3)," (",round(sd(beta_est_Raps),3),")")
mean(beta_cover_Raps)
sd(beta_est_Raps)









beta_est = rep(0,n.rep)
beta_cover = rep(0,n.rep)
beta_se = rep(0,n.rep)


beta_est_Raps = rep(0,n.rep)
beta_cover_Raps = rep(0,n.rep)

beta_est_IVW = rep(0,n.rep)
beta_cover_IVW = rep(0,n.rep)

beta_est_raps = rep(0,n.rep)
beta_cover_raps = rep(0,n.rep)

#balanced pleotropic effect
library(MASS)
for(k in 1:n.rep){
  if(k%%100){print(k)}
  R = diag(n.iv)
  #R[upper.tri(R)] = R[lower.tri(R)] = 0.5
  #mean_alpha = c(eff_alpha,rep(eff_alpha*0.5,n.iv-1))
  mean_alpha = rep(eff_alpha,n.iv)
  alpha = as.vector(rmvnorm(1,mu = mean_alpha,sigma = R*sigma_alpha))
  #alpha[1] = rnorm(1,mean_alpha,sd = sqrt(sigma_alpha))
  var_alpha = rep(sigma_alpha,n.iv)
  alpha_temp = rmvnorm(1,mu = mean_alpha,sigma = R*sigma_alpha)
  #alpha_temp[1] = rnorm(1,mean_alpha,sd = sqrt(sigma_alpha))
  error_y = rmvnorm(1,mu = rep(0,n.iv),sigma = R*sigma_Gamma)
  error_pleo = rnorm(n.iv,mean = 0, sd = sqrt(sigma_tau))
  Gamma = as.vector(alpha_temp*beta + error_pleo + error_y)
  var_Gamma = rep(sigma_Gamma+sigma_tau,n.iv)
  #R = diag(n.iv)
  MR_result <- MRWeight(Gamma,
                        var_Gamma,
                        alpha,
                        var_alpha,R)
  beta_est[k] = MR_result[[1]]
  beta_cover[k] = ifelse(MR_result[[2]]<=beta&MR_result[[3]]>=beta,1,0)
  beta_se[k] = MR_result[[4]]
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                      alpha,var_alpha)
  beta_est_IVW[k] = IVW_c_temp[1]
  beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
                               IVW_c_temp[3]>=beta,1,0)
  raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
                                           beta.outcome = Gamma,
                                           se.exposure = sqrt(var_alpha),
                                           se.outcome = sqrt(var_Gamma)),
                         diagnostics = F)
  
  beta_est_Raps[k] = raps_result$beta.hat
  beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                                raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
  # IVW_c_temp <- IVW_c(Gamma,var_Gamma,
  #                     alpha,var_alpha)
  # beta_est_IVW[k] = IVW_c_temp[1]
  # beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
  #                           IVW_c_temp[3]>=beta,1,0)
  # 
  # 
  raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
                                           beta.outcome = Gamma,
                                           se.exposure = sqrt(var_alpha),
                                           se.outcome = sqrt(var_Gamma)),
                         diagnostics = F)

  beta_est_Raps[k] = raps_result$beta.hat
  beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                             raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
  
}

paste0(round(mean(beta_est),3)," (",round(sd(beta_est),3),")")
mean(beta_cover)
mean(beta_se)
sd(beta_est)

paste0(round(mean(beta_est_IVW),3)," (",round(sd(beta_est_IVW),3),")")
mean(beta_cover_IVW)


paste0(round(mean(beta_est_Raps),3)," (",round(sd(beta_est_Raps),3),")")
mean(beta_cover_Raps)




beta_est = rep(0,n.rep)
beta_cover = rep(0,n.rep)
beta_se = rep(0,n.rep)


beta_est_Raps = rep(0,n.rep)
beta_cover_Raps = rep(0,n.rep)

beta_est_IVW = rep(0,n.rep)
beta_cover_IVW = rep(0,n.rep)

#inbalanced pleotropic effect
library(MASS)
for(k in 1:n.rep){
  if(k%%100){print(k)}
  R = diag(n.iv)
  #R[upper.tri(R)] = R[lower.tri(R)] = 0.5
  #mean_alpha = c(eff_alpha,rep(eff_alpha*0.5,n.iv-1))
  mean_alpha = rep(eff_alpha,n.iv)
  alpha = as.vector(rmvnorm(1,mu = mean_alpha,sigma = R*sigma_alpha))
  #alpha[1] = rnorm(1,mean_alpha,sd = sqrt(sigma_alpha))
  var_alpha = rep(sigma_alpha,n.iv)
  alpha_temp = rmvnorm(1,mu = mean_alpha,sigma = R*sigma_alpha)
  #alpha_temp[1] = rnorm(1,mean_alpha,sd = sqrt(sigma_alpha))
  error_y = rmvnorm(1,mu = rep(0,n.iv),sigma = R*sigma_Gamma)
  error_pleo = c(rep(0,n.iv/2),
    rnorm(n.iv/2,mean = 0, sd = sqrt(sigma_tau)))
  
  Gamma = as.vector(alpha_temp*beta + error_pleo + error_y)
  var_Gamma = c(rep(sigma_Gamma,n.iv/2),
    rep(sigma_Gamma+sigma_tau,n.iv/2))
  #R = diag(n.iv)
  MR_result <- MRWeight(Gamma,
                        var_Gamma,
                        alpha,
                        var_alpha,R)
  beta_est[k] = MR_result[[1]]
  beta_cover[k] = ifelse(MR_result[[2]]<=beta&MR_result[[3]]>=beta,1,0)
  beta_se[k] = MR_result[[4]]
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                      alpha,var_alpha)
  beta_est_IVW[k] = IVW_c_temp[1]
  beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
                               IVW_c_temp[3]>=beta,1,0)
  raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
                                           beta.outcome = Gamma,
                                           se.exposure = sqrt(var_alpha),
                                           se.outcome = sqrt(var_Gamma)),
                         diagnostics = F)
  
  beta_est_Raps[k] = raps_result$beta.hat
  beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                                raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
  
  # IVW_c_temp <- IVW_c(Gamma,var_Gamma,
  #                     alpha,var_alpha)
  # beta_est_IVW[k] = IVW_c_temp[1]
  # beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
  #                           IVW_c_temp[3]>=beta,1,0)
  # 
  # 
  # raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
  #                                          beta.outcome = Gamma,
  #                                          se.exposure = sqrt(var_alpha),
  #                                          se.outcome = sqrt(var_Gamma)),
  #                        diagnostics = F)
  # 
  # beta_est_Raps[k] = raps_result$beta.hat
  # beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
  #                            raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
  
}

paste0(round(mean(beta_est),3)," (",round(sd(beta_est),3),")")
mean(beta_cover)
mean(beta_se)
sd(beta_est)

paste0(round(mean(beta_est_IVW),3)," (",round(sd(beta_est_IVW),3),")")
mean(beta_cover_IVW)

paste0(round(mean(beta_est_Raps),3)," (",round(sd(beta_est_Raps),3),")")
mean(beta_cover_Raps)



ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

ar1_cor(n.iv,0.5)
library(susieR)
error_sigma = matrix(c(1,0.3,0.3,1),2,2)
#correlated SNPs
library(MASS)
n.rep = 1000
beta_est = rep(0,n.rep)
beta_cover = rep(0,n.rep)
beta_se = rep(0,n.rep)


beta_est_Raps = rep(0,n.rep)
beta_cover_Raps = rep(0,n.rep)

beta_est_IVW = rep(0,n.rep)
beta_cover_IVW = rep(0,n.rep)

for(k in 1:n.rep){
  if(k%%100){print(k)}
  R =ar1_cor(n.iv,0.5)
  G1 = rmvnorm(N,mu = rep(0,n.iv),R)
  error = rnorm(N,0,sd = 1)
  mean_alpha = rep(eff_alpha,n.iv)
  M1 = G1%*%mean_alpha+error
  sumstats <- univariate_regression(G1, M1)
  alpha = sumstats$betahat
  var_alpha = sumstats$sebetahat^2
  z = alpha/sqrt(var_alpha)
  p = 2*pnorm(-abs(z),lower.tail = T)
  
  G2 = rmvnorm(N,mu = rep(0,n.iv),R)
  
  error = rmvnorm(N,mu = rep(0,2),sigma = error_sigma)
  M2 = G2%*%mean_alpha+error[,1]
  Y2= M2*beta+error[,2]
  sumstats <- univariate_regression(G2, Y2)
  Gamma = sumstats$betahat
  var_Gamma = sumstats$sebetahat^2
  
  
  #R = diag(n.iv)
  MR_result <- MRWeight(Gamma,
                        var_Gamma,
                        alpha,
                        var_alpha,R)
  beta_est[k] = MR_result[[1]]
  beta_cover[k] = ifelse(MR_result[[2]]<=beta&MR_result[[3]]>=beta,1,0)
  beta_se[k] = MR_result[[4]]
  id.select = seq(5,n.iv,5)
  Gamma = Gamma[id.select]
  var_Gamma = var_Gamma[id.select]
  alpha = alpha[id.select]
  var_alpha = var_alpha[id.select]
  
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                      alpha,var_alpha)
  beta_est_IVW[k] = IVW_c_temp[1]
  beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
                               IVW_c_temp[3]>=beta,1,0)
  
  raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
                                           beta.outcome = Gamma,
                                           se.exposure = sqrt(var_alpha),
                                           se.outcome = sqrt(var_Gamma)),
                         diagnostics = F)

  beta_est_Raps[k] = raps_result$beta.hat
  beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                             raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
  # IVW_c_temp <- IVW_c(Gamma,var_Gamma,
  #                     alpha,var_alpha)
  # beta_est_IVW[k] = IVW_c_temp[1]
  # beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
  #                           IVW_c_temp[3]>=beta,1,0)
  # 
  # 
  # raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
  #                                          beta.outcome = Gamma,
  #                                          se.exposure = sqrt(var_alpha),
  #                                          se.outcome = sqrt(var_Gamma)),
  #                        diagnostics = F)
  # 
  # beta_est_Raps[k] = raps_result$beta.hat
  # beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
  #                            raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
  
}


paste0(round(mean(beta_est),3)," (",round(sd(beta_est),3),")")
mean(beta_cover)
mean(beta_se)
sd(beta_est)

paste0(round(mean(beta_est_IVW),3)," (",round(sd(beta_est_IVW),3),")")
mean(beta_cover_IVW)


paste0(round(mean(beta_est_Raps),3)," (",round(sd(beta_est_Raps),3),")")
mean(beta_cover_Raps)







library(susieR)
error_sigma = matrix(c(1,0.3,0.3,1),2,2)
#correlated SNPs
library(MASS)
n.rep = 1000
beta_est = rep(0,n.rep)
beta_cover = rep(0,n.rep)
beta_se = rep(0,n.rep)


beta_est_Raps = rep(0,n.rep)
beta_cover_Raps = rep(0,n.rep)

beta_est_IVW = rep(0,n.rep)
beta_cover_IVW = rep(0,n.rep)
p.causal = 0.2
for(k in 1:n.rep){
  if(k%%100){print(k)}
  #R =ar1_cor(n.iv,0.5)
  R =diag(n.iv)
  G1 = rmvnorm(N,mu = rep(0,n.iv),R)
  error = rnorm(N,0,sd = 1)
  mean_alpha = c(rep(eff_alpha,n.iv*p.causal),rep(0,n.iv*(1-p.causal)))
  M1 = G1%*%mean_alpha+error
  sumstats <- univariate_regression(G1, M1)
  alpha = sumstats$betahat
  var_alpha = sumstats$sebetahat^2
  z = alpha/sqrt(var_alpha)
  p = 2*pnorm(-abs(z),lower.tail = T)
  
  G2 = rmvnorm(N,mu = rep(0,n.iv),R)
  
  error = rmvnorm(N,mu = rep(0,2),sigma = error_sigma)
  M2 = G2%*%mean_alpha+error[,1]
  Y2= M2*beta+error[,2]
  sumstats <- univariate_regression(G2, Y2)
  Gamma = sumstats$betahat
  var_Gamma = sumstats$sebetahat^2
  
  
  #R = diag(n.iv)
  MR_result <- MRWeight(Gamma,
                        var_Gamma,
                        alpha,
                        var_alpha,R)
  beta_est[k] = MR_result[[1]]
  beta_cover[k] = ifelse(MR_result[[2]]<=beta&MR_result[[3]]>=beta,1,0)
  beta_se[k] = MR_result[[4]]
  #id.select = seq(5,n.iv,5)
  id.select = c(1:n.iv)
  Gamma = Gamma[id.select]
  var_Gamma = var_Gamma[id.select]
  alpha = alpha[id.select]
  var_alpha = var_alpha[id.select]
  
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                      alpha,var_alpha)
  beta_est_IVW[k] = IVW_c_temp[1]
  beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
                               IVW_c_temp[3]>=beta,1,0)
  
  raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
                                           beta.outcome = Gamma,
                                           se.exposure = sqrt(var_alpha),
                                           se.outcome = sqrt(var_Gamma)),
                         diagnostics = F)
  
  beta_est_Raps[k] = raps_result$beta.hat
  beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                                raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
  # IVW_c_temp <- IVW_c(Gamma,var_Gamma,
  #                     alpha,var_alpha)
  # beta_est_IVW[k] = IVW_c_temp[1]
  # beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
  #                           IVW_c_temp[3]>=beta,1,0)
  # 
  # 
  # raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
  #                                          beta.outcome = Gamma,
  #                                          se.exposure = sqrt(var_alpha),
  #                                          se.outcome = sqrt(var_Gamma)),
  #                        diagnostics = F)
  # 
  # beta_est_Raps[k] = raps_result$beta.hat
  # beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
  #                            raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
  
}


paste0(round(mean(beta_est),3)," (",round(sd(beta_est),3),")")
mean(beta_cover)
mean(beta_se)
sd(beta_est)

paste0(round(mean(beta_est_IVW),3)," (",round(sd(beta_est_IVW),3),")")
mean(beta_cover_IVW)


paste0(round(mean(beta_est_Raps),3)," (",round(sd(beta_est_Raps),3),")")
mean(beta_cover_Raps)

























n.rep = 3000

n.iv = 30

sigma_tau = 0.01^2

beta_est = rep(0,n.rep)
beta_cover = rep(0,n.rep)
beta_se = rep(0,n.rep)


beta_est_Raps = rep(0,n.rep)
beta_cover_Raps = rep(0,n.rep)

beta_est_IVW = rep(0,n.rep)
beta_cover_IVW = rep(0,n.rep)

beta_est_raps = rep(0,n.rep)
beta_cover_raps = rep(0,n.rep)

p.causal = 0.2
#independent IVs, no pleotropic
for(k in 1:n.rep){
  if(k%%100){print(k)}
  R = diag(n.iv)
  #R[upper.tri(R)] = R[lower.tri(R)] = 0.5
  #mean_alpha = c(eff_alpha,rep(eff_alpha*0.5,n.iv-1))
  mean_alpha = c(rep(eff_alpha,n.iv*p.causal),rep(0,n.iv*(1-p.causal)))
  
  
  
  alpha = as.vector(rmvnorm(1,mu = mean_alpha,sigma = R*sigma_alpha))
  #alpha[1] = rnorm(1,mean_alpha,sd = sqrt(sigma_alpha))
  var_alpha = rep(sigma_alpha,n.iv)
  alpha_temp = rmvnorm(1,mu = mean_alpha,sigma = R*sigma_alpha)
  #alpha_temp[1] = rnorm(1,mean_alpha,sd = sqrt(sigma_alpha))
  error_y = rmvnorm(1,mu = rep(0,n.iv),sigma = R*sigma_Gamma)
  Gamma = as.vector(alpha_temp*beta  + error_y)
  var_Gamma = rep(sigma_Gamma,n.iv)
  #R = diag(n.iv)
  MR_result <- MRWeight(Gamma,
                        var_Gamma,
                        alpha,
                        var_alpha,R)
  beta_est[k] = MR_result[[1]]
  beta_cover[k] = ifelse(MR_result[[2]]<=beta&MR_result[[3]]>=beta,1,0)
  beta_se[k] = MR_result[[4]]
  
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                      alpha,var_alpha)
  beta_est_IVW[k] = IVW_c_temp[1]
  beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
                               IVW_c_temp[3]>=beta,1,0)
  
  # raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
  #                                          beta.outcome = Gamma,
  #                                          se.exposure = sqrt(var_alpha),
  #                                          se.outcome = sqrt(var_Gamma)),
  #                        diagnostics = F)
  # 
  # beta_est_raps[k] = raps_est = raps_result$beta.hat
  # raps_se = raps_result$beta.se
  # beta_cover_raps[k] = ifelse(raps_est-1.96*raps_se<=beta&
  #                               raps_est+1.96*raps_se>=beta,1,0)
  
  # 
  # 
  raps_result <- mr.raps(data = data.frame(beta.exposure = alpha,
                                           beta.outcome = Gamma,
                                           se.exposure = sqrt(var_alpha),
                                           se.outcome = sqrt(var_Gamma)),
                         diagnostics = F,
                         over.dispersion = F,
                         loss.function = "l2",
                         se.method = "sandwich")
  
  beta_est_Raps[k] = raps_result$beta.hat
  beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                                raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
  
}
# cov(beta_est,beta_cover)
# 
# 0.5*(sigma_Gamma+0.15^2*sigma_alpha)


paste0(round(mean(beta_est),3)," (",round(sd(beta_est),3),")")
mean(beta_cover)
mean(beta_se)
sd(beta_est)

paste0(round(mean(beta_est_IVW),3)," (",round(sd(beta_est_IVW),3),")")
mean(beta_cover_IVW)
paste0(round(mean(beta_est_Raps),3)," (",round(sd(beta_est_Raps),3),")")
mean(beta_cover_Raps)
sd(beta_est_Raps)




