#Get the coverage probability for single ratio distribution
#Three different method were tested
#1. Directly use ratio estimate
#2. Directly use ratio estimate with second term
#Ratio estimate
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])
i3 = as.numeric(args[[3]])
i4 = as.numeric(args[[4]])
setwd("/data/zhangh24/MR_MA/")
Regression = function(Y,M,G,G2){
  n = length(Y)
  model1 = lm(Y~G-1)
  coef_temp = coef(summary(model1))
  Gamma = coef_temp[1]
  var_Gamma = coef_temp[2]^2
  model2 = lm(M~G2-1)
  coef_temp2 = coef(summary(model2))
  gamma = coef_temp2[1]
  var_gamma = coef_temp2[2]^2
  return(c(Gamma,var_Gamma,
           gamma,var_gamma))
  
  
}
#standard delta method


beta_vec <- c(0,0.2,0.5)
n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.00,0.01,0.03,0.05)

times = 1000
n <- n_vec[i1]
MAF =0.25
Gamma_est <- rep(0,times)
Gamma_var <- rep(0,times)
gamma_est <- rep(0,times)
gamma_var <- rep(0,times)

ratio_var_standard <- rep(0,times)
G_ori = matrix(rbinom(n*5,1,MAF),n,1)
G = apply(G_ori,2,scale)
G_ori2 = matrix(rbinom(n*5,1,MAF),n,1)
G2 = apply(G_ori2,2,scale)
set.seed(i3)
for(i in 1:times){
  print(i)
  beta_M = beta_vec[i4]
  alpha_G = alpha_vec[i2]
  sigma_y = 1
  sigma_m = 1
  U = rnorm(n,sd = 1)
  alpha_U <- 0.1
  beta_U <- 0.1
  M = G%*%alpha_G + alpha_U*U+rnorm(n,sd = sqrt(sigma_m))
  Y = beta_M*M +beta_U*U+rnorm(n,sd = sqrt(sigma_y))
  U2 = rnorm(n,sd = 1)
  M = G2%*%alpha_G + alpha_U*U2+rnorm(n,sd = sqrt(sigma_m))
  est <- Regression(Y,M,G,G2)
  
  Gamma_est[i] <- est[1]
  Gamma_var[i] <- est[2]
  gamma_est[i] <- est[3]
  gamma_var[i] <- est[4]
  
}
ratio_est = Gamma_est/gamma_est

#just using the leading term
ratio_var_standard = Gamma_var/gamma_est^2

z_est <- ratio_est/sqrt(var_ratio_standard)
p_est <- 2*pnorm(-abs(z_est))
ci_low_ratio <- ratio_est-sqrt(ratio_var_standard)*1.96
ci_high_ratio <- ratio_est+sqrt(ratio_var_standard)*1.96
cover_ratio <- ifelse(beta_M>=ci_low_ratio&
                        beta_M<=ci_high_ratio,1,0)
#add the second term
ratio_var_alpha = Gamma_var/gamma_est^2+ gamma_var*Gamma_est^2/gamma_est^4
ci_low_alpha <- ratio_est-sqrt(ratio_var_alpha)*1.96
ci_high_alpha <- ratio_est+sqrt(ratio_var_alpha)*1.96
cover_alpha <- ifelse(beta_M>=ci_low_alpha&
                        beta_M<=ci_high_alpha,1,0)

result <- list(Gamma_est,
               Gamma_var,
               gamma_est,
               gamma_var,
               ratio_est,
               ratio_var_standard,
               ratio_var_alpha,
               cover_true,
               cover_ratio,
               cover_alpha,
               ci_low_ratio,
               ci_high_ratio,
               ci_low_alpha,
               ci_high_alpha)

save(result,file = paste0("./result/simulation/ratio_estimate/ratio_estimate_delta_",i1,"_",i2,"_",i3,"_",i4,".Rdata"))




