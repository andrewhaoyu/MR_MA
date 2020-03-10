#Get the coverage probability for single ratio distribution
#Three different method were tested
#1. Directly use ratio estimate
#2. Use Zhonghua's Z statistics method (epi)
#3. Use exact distribution method (exact)
#method 1 and 2 are implemented in Ratio estimate
#method 3 is implemented in RatioExact
#Ratio estimate
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])
i3 = as.numeric(args[[3]])
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

Ratio = function(Gamma,var_Gamma,gamma,var_gamma,n){
  ratio_est = Gamma/gamma
  var_ratio = var_Gamma/gamma^2+var_gamma*Gamma^2/gamma^4
  n.simu <- 1000000
  z_Gamma <- rnorm(n.simu)
  z_gamma <- rnorm(n.simu,mean = gamma*sqrt(n),sd = 1)
  true_distribution <- z_Gamma/sqrt(1+z_Gamma^2/z_gamma^2)
  q_result <- quantile(true_distribution,c(0.025,0.975))
  ci_low <- ratio_est-q_result[2]*sqrt(var_ratio)
  ci_high <- ratio_est-q_result[1]*sqrt(var_ratio)
  cover = ifelse(beta_M>=ci_low&
                   beta_M<=ci_high,1,0)
  return(c(ratio_est,var_ratio,cover,ci_low,ci_high))  
}
#new ratio formula
# RatioExact = function(Gamma,var_Gamma,gamma,var_gamma,n){
#   ratio_est = Gamma/gamma
#   n.simu <- 1000000
#   z_Gamma <- rnorm(n.simu,mean =0,sd =sqrt(var_Gamma))
#   z_gamma <- rnorm(n.simu,mean = gamma,sd = sqrt(var_gamma))
#   true_distribution <- z_Gamma/z_gamma
#   q_result <- quantile(true_distribution,c(0.025,0.975))
#   cover = ifelse(ratio_est>=q_result[1]&
#                    ratio_est<=q_result[2],1,0)
#   ci_low <- ratio_est-q_result[2]
#   ci_high <- ratio_est-q_result[1]
#   return(c(cover,ci_low,ci_high))
# }

RatioExact = function(Gamma,var_Gamma,gamma,var_gamma,n){
  # i = idx[1]
  # 
  #  Gamma <-  Gamma_est[i]
  #  var_Gamma = Gamma_var[i]
  #  gamma = gamma_est[i]
  # var_gamma = gamma_var[i]
  ratio_est = Gamma/gamma
  n.simu <- 1000000
   # z_Gamma <- rnorm(n.simu,mean =0,sd =sqrt(var_Gamma))
   # z_gamma <- rnorm(n.simu,mean = gamma,sd = sqrt(var_gamma))
  #z_Gamma <- rnorm(n.simu,mean =0,sd =sqrt(var_Gamma))
  #z_gamma <- rnorm(n.simu,mean = gamma,sd = sqrt(var_gamma))
   #z_gamma <- rnorm(n.simu,mean = sqrt(n)*alpha_G,sd = sqrt((n-1)*var_gamma))
   z_Gamma <- rnorm(n.simu,mean =0,sd =sqrt(var_Gamma))
   z_gamma <- rnorm(n.simu,mean = gamma,sd = sqrt(var_gamma))
  true_distribution <- z_Gamma/z_gamma
  q_result <- quantile(true_distribution,c(0.025,0.975))
  cover = ifelse(ratio_est>=q_result[1]&
                   ratio_est<=q_result[2],1,0)
  #the q_result should always be bounded by quachy distribution
  #q_result[1] = ifelse(q_result[1]<=qcauchy(0.025),qcauchy(0.025),q_result[1])
  #q_result[2] = ifelse(q_result[2]>=qcauchy(0.975),qcauchy(0.975),q_result[2])
  ci_low <- ratio_est-q_result[2]
  ci_high <- ratio_est-q_result[1]
 
  return(c(cover,ci_low,ci_high))
}
  
# idx.temp <- which(cover_exact!=
#                     cover_true_exact)
# Gamma = Gamma_est[idx.temp[1]]
# var_Gamma = Gamma_var[idx.temp[1]]
# gamma = gamma_est[idx.temp[1]]
# var_gamma = gamma_var[idx.temp[1]]
# RatioExact = function(Gamma,var_Gamma,gamma,var_gamma,n){
#   ratio_est = Gamma/gamma
#   n.simu <- 1000000
#   z_Gamma <- rt(n.simu,df =n-1)*sqrt((n-1)*var_Gamma)
#   z_gamma <- rt(n.simu,df = n-1,ncp = sqrt(n)*gamma/(sqrt((n-1)*var_gamma)))*sqrt((n-1)*var_gamma)
#   true_distribution <- z_Gamma/z_gamma
#   q_result <- quantile(true_distribution,c(0.025,0.975))
#   cover = ifelse(ratio_est>=q_result[1]&
#                    ratio_est<=q_result[2],1,0)
#   ci_low <- ratio_est-q_result[2]
#   ci_high <- ratio_est-q_result[1]
#   return(c(cover,ci_low,ci_high))
# }






n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.00,0.01,0.03,0.05)

times = 1000
n <- n_vec[i1]
MAF =0.25
Gamma_est <- rep(0,times)
Gamma_var <- rep(0,times)
gamma_est <- rep(0,times)
gamma_var <- rep(0,times)
ratio_est <- rep(0,times)
ratio_var <- rep(0,times)
ratio_cover <- rep(0,times)
cover_epi <- rep(0,times)
cover_exact <- rep(0,times)
ci_low_epi <- rep(0,times)
ci_high_epi <- rep(0,times)
ci_low_exact <- rep(0,times)
ci_high_exact <- rep(0,times)
G_ori = matrix(rbinom(n*5,1,MAF),n,1)
G = apply(G_ori,2,scale)
G_ori2 = matrix(rbinom(n*5,1,MAF),n,1)
G2 = apply(G_ori2,2,scale)
set.seed(i3)
for(i in 1:times){
  print(i)
  beta_M = 0
  alpha_G = alpha_vec[i2]
  sigma_y = 1
  sigma_m = 1
  U = rnorm(n,sd = 1)
  alpha_U <- 0.1
  beta_U <- 0.1
  M = G%*%alpha_G + rnorm(n,sd = sqrt(sigma_m))
  Y = beta_M*M +rnorm(n,sd = sqrt(sigma_y))
  U2 = rnorm(n,sd = 1)
  M = G2%*%alpha_G + rnorm(n,sd = sqrt(sigma_m))
  est <- Regression(Y,M,G,G2)
  
  Gamma_est[i] <- est[1]
  Gamma_var[i] <- est[2]
  gamma_est[i] <- est[3]
  gamma_var[i] <- est[4]
  ratio_temp <- Ratio(est[1],est[2],est[3],est[4],n)
  ratio_est[i] <- ratio_temp[1]
  ratio_var[i] <- ratio_temp[2]
  cover_epi[i] <- ratio_temp[3]
  ci_low_epi[i] <- ratio_temp[4]
  ci_high_epi[i] <- ratio_temp[5]
  ratio_exact_temp <- RatioExact(est[1],est[2],est[3],est[4],n)
  cover_exact[i] <- ratio_exact_temp[1]
  ci_low_exact[i] <- ratio_exact_temp[2]
  ci_high_exact[i] <- ratio_exact_temp[3]
  
}


idx <- which(cover_exact==0)



n.simu <- 1000000
z_est <- ratio_est/sqrt(ratio_var)
p_est <- 2*pnorm(-abs(z_est))
ci_low_ratio <- ratio_est-sqrt(ratio_var)*1.96
ci_high_ratio <- ratio_est+sqrt(ratio_var)*1.96
cover_ratio <- ifelse(beta_M>=ci_low_ratio&
                        beta_M<=ci_high_ratio,1,0)


standard_norm <- rnorm(times)
z_Gamma <- rnorm(n.simu)
z_gamma <- rnorm(n.simu,mean = alpha_G*sqrt(n),sd = 1)
true_distribution <- z_Gamma/sqrt(1+z_Gamma^2/z_gamma^2)
q_result <- quantile(true_distribution,c(0.025,0.975))


cover_true = ifelse(z_est>=q_result[1]&
                 z_est<=q_result[2],1,0)
#cover_true <- sum(cover_vec)/length(cover_vec)

z_Gamma <- rnorm(n.simu,mean = 0, sd = sqrt(sigma_y/n))
z_gamma <- rnorm(n.simu,mean = alpha_G,sd = sqrt(sigma_m/n))

true_distribution <- z_Gamma/z_gamma
q_result <- quantile(true_distribution,c(0.025,0.975))
q_result[1] = ifelse(q_result[1]<=qcauchy(0.025),qcauchy(0.025),q_result[1])
q_result[2] = ifelse(q_result[2]>=qcauchy(0.975),qcauchy(0.975),q_result[2])

print(q_result)
cover_true_exact = ifelse(ratio_est>=q_result[1]&
                           ratio_est<=q_result[2],1,0)




#cover_epi <- sum(cover_epi)/length(cover_epi)
#cover_exact <- sum(cover_exact)/length(cover_exact)
#cover_true_exact <- sum(cover_vec_exact)/length(cover_vec_exact)

result <- list(Gamma_est,
               Gamma_var,
               gamma_est,
               gamma_var,
               ratio_est,
               ratio_var,
               cover_ratio,
               cover_true,
               cover_epi,
               cover_exact,
               cover_true_exact,
               ci_low_ratio,
               ci_high_ratio,
               ci_low_epi,
               ci_high_epi,
               ci_low_exact,
               ci_high_exact)

save(result,file = paste0("./result/simulation/ratio_estimate/ratio_estimate_",i1,"_",i2,"_",i3,".Rdata"))




