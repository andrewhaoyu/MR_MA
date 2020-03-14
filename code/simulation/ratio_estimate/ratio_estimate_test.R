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
  
  # if((gamma)<=1.96*sqrt(var_gamma)&(gamma)>=-1.96*sqrt(var_gamma)){
  #   plug_mean = 0
  # }else{
  #   plug_mean = gamma*sqrt(n)
  # }
  z_gamma <- rnorm(n.simu,mean = gamma*sqrt(n),sd = 1)
  true_distribution <- z_Gamma/sqrt(1+z_Gamma^2/z_gamma^2)
 T_ratio <- ratio_est/sqrt(var_ratio)
  p = 2*sum(true_distribution<=-abs(T_ratio))/n.simu
  
 
  return(c(ratio_est,var_ratio,p))  
}


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
ratio_p1 = rep(0,times)
ratio_var_standard <- rep(0,times)
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
  #estimate the probability if gamma is plugged in
  ratio_p1[i] <- ratio_temp[3]
  
}
T_ratio = ratio_est/sqrt(ratio_var)
#estimate the probability if gamma is 0
n.simu <- 1000000
z_Gamma <- rnorm(n.simu)
z_gamma <- rnorm(n.simu,mean = 0,sd = 1)
true_distribution <- z_Gamma/sqrt(1+z_Gamma^2/z_gamma^2)

for(i in 1:times){
  ratio_p2[i] = 2*sum(true_distribution<=-abs(T_ratio[i]))/n.simu
  
}

#estimate the proportion of gamma to be nonzero
z_gamma = gamma_est/sqrt(gamma_var)
#number of null
p0=sum(z_gamma>=-1.96&z_gamma<=1.96)/length(z_gamma)
p_weight = p0*ratio_p2+(1-p0)*ratio_p1
sum(p_weight<=0.05)/length(p_weight)
sum(ratio_p1<=0.05)/length(ratio_p1)


var_ratio_standard = Gamma_var/gamma_est^2
z_est <- ratio_est/sqrt(var_ratio_standard)
p_delta <- 2*pnorm(-abs(z_est))
sum(p_delta<=0.05)/length(p_delta)
result <- list(Gamma_est,
               Gamma_var,
               gamma_est,
               gamma_var,
               ratio_est,
               ratio_var,
               ratio_p1,
               ratio_p2,
               p_weight,
               p_delta)

save(result,file = paste0("./result/simulation/ratio_estimate/p_estimate_",i1,"_",i2,"_",i3,".Rdata"))




