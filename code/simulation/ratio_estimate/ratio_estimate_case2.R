#Get the coverage probability for single ratio distribution
#Ratio estimate
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])
i3 = as.numeric(args[[3]])
setwd("/data/zhangh24/MR_MA/")
Regression = function(Y,M,G){
  n = length(Y)
  model1 = lm(Y~G-1)
  coef_temp = coef(summary(model1))
  Gamma = coef_temp[1]
  var_Gamma = coef_temp[2]^2
  model2 = lm(M~G-1)
  coef_temp2 = coef(summary(model2))
  gamma = coef_temp2[1]
  var_gamma = coef_temp2[2]^2
  return(c(Gamma,var_Gamma,
           gamma,var_gamma))
  
  
}

Ratio = function(Gamma,var_Gamma,gamma,var_gamma,n){
  ratio_est = Gamma/gamma
  var_ratio = var_Gamma/gamma^2+var_gamma*Gamma^2/gamma^4
  n.simu <- 10000
  z_Gamma <- rnorm(n.simu)
  z_gamma <- rnorm(n.simu,mean = gamma*sqrt(n),sd = 1)
  z_est = ratio_est/sqrt(var_ratio)
  p = pchisq(4*z_est^2,df=1,lower.tail=F)
  cover = ifelse(p>=0.05,1,0)
  
  ci_low <- -sqrt(qchisq(p=0.95,df=1)/4)*sqrt(var_ratio)+ratio_est
  ci_high <- sqrt(qchisq(p=0.95,df=1)/4)*sqrt(var_ratio)+ratio_est
  
  return(c(ratio_est,var_ratio,cover,
           ci_low,ci_high))  
}

set.seed(i3)
n_vec <- c(15000,75000,150000)
beta_vec <- c(0)

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
ci_low_epi <- rep(0,times)
ci_high_epi <- rep(0,times)
G_ori = matrix(rbinom(n*5,1,MAF),n,1)
G = apply(G_ori,2,scale)
for(i in 1:times){
  #print(i)
  beta_M = 0
  beta_G = beta_vec[i2]
  sigma_y = 1
  sigma_m = 1
  M = G%*%beta_G+rnorm(n,sd = sqrt(sigma_m))
  Y = rnorm(n,sd = sqrt(sigma_y))
  est <- Regression(Y,M,G)
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
}

z_est <- ratio_est/sqrt(ratio_var)
standard_norm <- rnorm(times)
z_Gamma <- rnorm(times)
z_gamma <- rnorm(times,mean = beta_G*sqrt(n),sd = 1)

true_distribution <- z_Gamma/sqrt(1+z_Gamma^2/z_gamma^2)
p_est <- 2*pnorm(-abs(z_est))
cover_ratio <- 1-sum(p_est<=0.05)/times
ci_low_ratio <- ratio_est-sqrt(ratio_var)*1.96
ci_high_ratio <- ratio_est+sqrt(ratio_var)*1.96

q_result <- quantile(true_distribution,c(0.025,0.975))
cover_vec = ifelse(z_est>=q_result[1]&
                     z_est<=q_result[2],1,0)

cover_true <- sum(cover_vec)/length(cover_vec)
cover_epi <- sum(cover_epi)/length(cover_epi)

result <- list(Gamma_est,
               Gamma_var,
               gamma_est,
               gamma_var,
               ratio_est,
               ratio_var,
               ratio_cover,
               cover_ratio,
               cover_true,
               cover_epi,
               ci_low_ratio,
               ci_high_ratio,
               ci_low_epi,
               ci_high_epi)
save(result,file = paste0("./result/simulation/ratio_estimate_case2",i1,"_",i2,"_",i3,".Rdata"))





# library(ggplot2)
# data <- data.frame(z_est,standard_norm,true_distribution)
# colnames(data) <- c("Ratio","Standard Normal","Derived_distribution")
# library(reshape2)
# data.m <- melt(data)
# ggplot(data.m,aes(value,colour=variable))+
#   geom_density()+
#   theme_Publication()



