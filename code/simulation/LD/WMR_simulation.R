#Goal: Simulate data with realistic LD information
#Y = M\beta + G\theta + U\beta_u + error_y
#M = G\alpha + U\alpha_U + error_m
#h2_y  = var(\beta G\alpha+G\theta) = 0.4
#h2_m = var(\alphaG) = 0.4
#var(error_m+U\alpha_U) = 0.6
#var(error_y+U\beta_U) = 0.6
#causal SNPs proportion for M: 0.1, 0.01
#overlapping between pleotripic and non pleotropic 1, 0.5, 0.75
n.snp = 100
beta = 0.15
N = 6000
cau.pro = 0.1
n.cau = n.snp*cau.pro
h2_m = 0.4
h2_y = 0.4
sigma_alpha = h2_m/n.cau
sigma_theta = (h2_y-beta^2)/n.cau
alpha_u = 0.547
sigma_error_m = 1-h2_m-alpha_u^2
beta_u = 0.547
sigma_error_y = 1-h2_y-beta_u^2

set.seed(123)
idx.cau_m = sample(c(1:n.snp),n.cau)
#plotropic settings
pleosnp.pro = 1
n.cau.overlap = pleosnp.pro*n.cau
n.cau.specific = n.cau - n.cau.overlap
#pleotrpic snps proportion the same as causal snps
idx.cau_pleo = c(sample(idx.cau_m,n.cau.overlap),
                 sample(setdiff(c(1:n.cau),idx.cau_m),n.cau-n.cau.overlap))

alpha = rnorm(n.cau,mean = 0,sd = sqrt(sigma_alpha))
theta = rnorm(n.cau,mean = 0, sd = sqrt(sigma_theta))
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
library(mr.raps)
library(Rfast)
library(MASS)
R =ar1_cor(n.snp,0.5)
G1 = rmvnorm(N,mu = rep(0,n.snp),R)
G2 = rmvnorm(N,mu = rep(0,n.snp),R)
U1 = rnorm(N)
U2 = rnorm(N)
G1.cau = G1[,idx.cau_m]
G2.cau = G2[,idx.cau_m]
G1.pleo = G1[,idx.cau_pleo]
G2.pleo = G2[,idx.cau_pleo]
#G1 to obtain sum data for Y
#G2 to obtain sum data for M
n.rep = 1000
beta_est = rep(0,n.rep)
beta_cover = rep(0,n.rep)
beta_se = rep(0,n.rep)


beta_est_Raps = rep(0,n.rep)
beta_cover_Raps = rep(0,n.rep)

beta_est_IVW = rep(0,n.rep)
beta_cover_IVW = rep(0,n.rep)

for(k in 1:n.rep){
  print(k)
  error_m = rnorm(N,sd = sqrt(sigma_error_m))
  error_y = rnorm(N,sd = sqrt(sigma_error_y))
  M1 = G1.cau%*%alpha+U1*alpha_u+error_m
  Y1 = M1%*%beta + G1.pleo%*%theta+U1*beta_u + error_y
  error_m = rnorm(N,sd = sqrt(sigma_error_m))
  M2 = G2.cau%*%alpha+U2*alpha_u+error_m
  
  library(susieR)
  sumstats <- univariate_regression(G1, M1)
  sumalpha = sumstats$betahat
  var_alpha = sumstats$sebetahat^2
  p_alpha = 2*pnorm(-abs(alpha/sqrt(var_alpha)),lower.tail = T)
  sumstats <- univariate_regression(G2, Y2)
  sumGamma = sumstats$betahat
  var_Gamma = sumstats$sebetahat^2
  #p_Gamma = 2*pnorm(-abs(Gamma/sqrt(var_Gamma)),lower.tail = T)
  Myclumping <- function(R,p){
    n.snp = ncol(R)
    
    
    #keep snps for clumpinp
    keep.ind = c(1:n.snp)
    #remove snps due to clumping
    remove.ind  = NULL
    #select snp ind
    select.ind = NULL
    temp = 1
    while(length(keep.ind)>0){
      # print(temp)
      p.temp = p[keep.ind]
      #select top snp
      top.ind = which.min(p.temp)
      select.ind = c(select.ind,keep.ind[top.ind])
      #print(keep.ind[top.ind])
      #tempory correlation
      R.temp = R[keep.ind[top.ind],]
      idx.remove = which(R.temp>=0.01)
      #take out
      remove.ind= c(remove.ind,idx.remove)
      keep.ind = setdiff(keep.ind,remove.ind)
      temp = temp+1
    }
    result = data.frame(select.ind,p[select.ind])
    return(result)
  }
  clump.snp = Myclumping(R,p_alpha)
  select.id = clump.snp[clump.snp$p.select.ind.<=5E-08,1]
  alpha_select =sumalpha[select.id]
  var_alpha_select = var_alpha[select.id]
  Gamma_select = sumGamma[select.id]
  var_Gamma_select = var_Gamma[select.id]
  
  IVW_c_temp <- IVW_c(Gamma_select,var_Gamma_select,
                      alpha_select,var_alpha_select)
  beta_est_IVW[k] = IVW_c_temp[1]
  beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
                               IVW_c_temp[3]>=beta,1,0)
  beta_est_IVW[k] = IVW_c_temp[1]
  beta_cover_IVW[k] = ifelse(IVW_c_temp[2]<=beta&
                               IVW_c_temp[3]>=beta,1,0)
  
  raps_result <- mr.raps(data = data.frame(beta.exposure = alpha_select,
                                           beta.outcome = Gamma_select,
                                           se.exposure = sqrt(var_alpha_select),
                                           se.outcome = sqrt(var_Gamma_select)),
                         diagnostics = F)
  beta_est_Raps[k] = raps_result$beta.hat
  beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                                raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  
  #R.select = R[select.id,select.id]
  MR_result <- MRWeight(Gamma = sumGamma,
                        var_Gamma = var_Gamma,
                        alpha = sumalpha,
                        var_alpha = var_alpha,
                        R = R)
  beta_est[k] = MR_result[[1]]
  beta_cover[k] = ifelse(MR_result[[2]]<=beta&MR_result[[3]]>=beta,1,0)
  beta_se[k] = MR_result[[4]]
  
}
paste0(round(mean(beta_est),3)," (",round(sd(beta_est),3),")")
mean(beta_cover)
mean(beta_se)
sd(beta_est)

paste0(round(mean(beta_est_IVW),3)," (",round(sd(beta_est_IVW),3),")")
mean(beta_cover_IVW)


paste0(round(mean(beta_est_Raps),3)," (",round(sd(beta_est_Raps),3),")")
mean(beta_cover_Raps)


# 
# MR_result <- MRWeight(Gamma = Gamma,
#                       var_Gamma = var_Gamma,
#                       gamma = gamma,
#                       var_gamma = var_gamma,
#                       R = R.select)


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




ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}

#ar1_cor(n.iv,0.5)

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




