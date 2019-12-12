#test the type one error and coverage prob of two-stage, IVW and IVW summary level statistics
#IVW estimate the variance using the meta-analysis variance
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])


setwd("/spin1/users/zhangh24/MR_MA/")
Meta = function(coef_vec,var_vec){
  meta_var = (sum(1/var_vec))^-1
  meta_coef = meta_var*sum(coef_vec/var_vec)
  return(c(meta_coef,meta_var))
}
MetaWeight = function(coef_vec,var_vec){
  meta_var = (sum(1/var_vec))^-1
  weight = meta_var/var_vec
  return(weight)
}




Project <- function(G){
  if(ncol(G)==1){
    G%*%t(G)/as.numeric(crossprod(G))
  }else{
    G%*%solve(t(G)%*%G)%*%t(G)  
  }
  
}
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
#IVW estimate using summary level statistics
IVW_s = function(Gamma,var_Gamma,gamma,var_gamma){
  p <- length(Gamma)
  raio_vec = rep(0,p)
  ratio_var_vec = rep(0,p)
  for(k in 1:p){
    raio_vec[k] = Gamma/gamma
    ratio_var_vec[k] =  var_Gamma/gamma^2+var_gamma*Gamma^2/gamma^4
  }
  
  Meta_result = Meta(raio_vec,ratio_var_vec)
  ratio_ivw =   Meta_result[1]
  ratio_ivw_var = Meta_result[2]
  coef_low = ratio_ivw-1.96*sqrt(ratio_ivw_var)
  coef_high = ratio_ivw+1.96*sqrt(ratio_ivw_var)
  cover = ifelse((beta_M>=coef_low&
                    beta_M<=coef_high),1,0)
  
  return(c(coef_est,cover,sigma_beta_est,
           coef_low,coef_high))
}



IVW_f = function(Y,M,G,beta_M){
  p <- length(Gamma)
  raio_vec = rep(0,p)
  ratio_var_vec = rep(0,p)
  for(k in 1:p){
    raio_vec[k] = Gamma/gamma
    ratio_var_vec[k] =  var_Gamma/gamma^2+var_gamma*Gamma^2/gamma^4
  }
  
  Meta_result = Meta(raio_vec,ratio_var_vec)
  ratio_ivw =   Meta_result[1]
  ratio_ivw_var = Meta_result[2]
  coef_low = ratio_ivw-1.96*sqrt(ratio_ivw_var)
  coef_high = ratio_ivw+1.96*sqrt(ratio_ivw_var)
  cover = ifelse((beta_M>=coef_low&
                    beta_M<=coef_high),1,0)
  
  Meta_Weight = Meta(coef_vec,var_vec)
  
  q_result <- quantile(true_distribution,c(0.025,0.975))
  z_est <- ratio_est/sqrt(var_ratio)
  cover = ifelse(z_est>=q_result[1]&
                   z_est<=q_result[2],1,0)
  ci_low <- q_result[1]*sqrt(var_ratio)+ratio_est
  ci_high <- q_result[2]*sqrt(var_ratio)+ratio_est
  
  coef_est =   Meta_result[1]
  sigma_beta_est = Meta_result[2]
  
  
  coef_low = coef_est-1.96*sqrt(sigma_beta_est)
  coef_high = coef_est+1.96*sqrt(sigma_beta_est)
  cover = ifelse((beta_M>=coef_low&
                    beta_M<=coef_high),1,0)
  
  return(c(coef_est,cover,sigma_beta_est))
}
GenerateTrueDistribution <- function( Gamma_vec,
var_Gamma_vec,
gamma_vec,
var_gamma_vec,
weight,
n,
var_vec){
 p <- length(Gamma_vec)
 n.simu <- 10000
 true_distribution <- rep(0,n.simu)
  for(k in 1:p){
    z_Gamma <- rnorm(n.simu)   
    z_gamma <- rnorm(n.simu,mean = gamma_vec[k]*sqrt(n),sd = 1)
    true_distribution_temp <- weight[k]*var_vec[k]*z_Gamma/sqrt(1+z_Gamma^2/z_gamma^2)
    true_distribution <- true_distribution+true_distribution_temp  
  }
 return(true_distribution)
}




n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.01,0.03,0.05)


set.seed(i1)
times = 100
n <- 15000
MAF =0.25
p <- 5
TwoStage_est = rep(0,times)
cover_TwoStage_est = rep(0,times)
sigma_TwoStage = rep(0,times)
sigma_y_TwoStage = rep(0,times)
IVW_est = rep(0,times)
cover_IVW_est = rep(0,times)
sigma_IVW = rep(0,times)
sigma_y_IVW = rep(0,times)
IVW_meta_est = rep(0,times)
cover_IVW_meta_est = rep(0,times)
sigma_IVW_meta = rep(0,times)
sigma_y_IVW_meta = rep(0,times)
IVWs_est = rep(0,times)
cover_IVWs_est = rep(0,times)
sigma_IVWs = rep(0,times)

#sigma_est = rep(0,times)
#sigma_beta_est = rep(0,times)
G_ori = matrix(rbinom(n*5,1,MAF),n,p)
G = apply(G_ori,2,scale)
for(i in 1:times){
  print(i)
  beta_M = 0.1
  beta_U = 0.1
  
  alpha_U = 0.01
  U = rnorm(n)
  # p = 5
  # MAF=0.25
  beta_G = rep(beta_G_vec[i2],p)
  sigma_y = 1
  sigma_m = 1
  M = G%*%beta_G+rnorm(n,sd = sqrt(sigma_m))
  Y = M*beta_M +rnorm(n,sd = sqrt(sigma_y))
  
  #M = G%*%beta_G+U*alpha_U+rnorm(n,sd = sqrt(sigma_m))
  #Y = M*beta_M + U*beta_U+rnorm(n,sd = sqrt(sigma_y))
  #M = G%*%beta_G+rnorm(n,sd = sqrt(sigma_m))
  #Y = M*beta_M + rnorm(n,sd = sqrt(sigma_y))
  TwoStage_result = TwoStage(Y,M,G,beta_M)  
  TwoStage_est[i] = TwoStage_result[1]
  cover_TwoStage_est[i] = TwoStage_result[2]
  sigma_TwoStage[i] = TwoStage_result[3]
  sigma_y_TwoStage[i] = TwoStage_result[4]
  IVW_result = IVW(Y,M,G,beta_M)  
  IVW_est[i] = IVW_result[1]
  cover_IVW_est[i] = IVW_result[2]
  sigma_IVW[i] = IVW_result[3]
  sigma_y_IVW[i] = IVW_result[4]
  IVW_meta_result = IVW_meta(Y,M,G,beta_M)  
  IVW_meta_est[i] = IVW_meta_result[1]
  cover_IVW_meta_est[i] = IVW_meta_result[2]
  sigma_IVW_meta[i] = IVW_meta_result[3]
  sigma_y_IVW_meta[i] = IVW_meta_result[4]
  IVWs_result = IVW_s(Y,M,G,beta_M)  
  IVWs_est[i] = IVWs_result[1]
  cover_IVWs_est[i] = IVWs_result[2]
  sigma_IVWs[i] = IVWs_result[3]
  
}

result1 = list(TwoStage_est,
               cover_TwoStage_est,
               sigma_TwoStage,
               sigma_y_TwoStage,
               IVW_est,
               cover_IVW_est,
               sigma_IVW,
               sigma_y_IVW,
               IVW_meta_est,
               cover_IVW_meta_est,
               sigma_IVW_meta,
               sigma_y_IVW_meta,
               IVWs_est,
               cover_IVWs_est,
               sigma_IVWs
)

mean(TwoStage_est)-beta_M
mean(IVW_est)-beta_M
mean(IVW_meta_est)-beta_M
mean(IVWs_est)-beta_M
mean(cover_TwoStage_est)
mean(cover_IVW_est)
mean(cover_IVW_meta_est,na.rm=T)
mean(cover_IVWs_est,na.rm = T)
var(TwoStage_est)
var(IVW_est)
var(IVW_meta_est)
var(IVWs_est)
mean(sigma_TwoStage)
mean(sigma_IVW)
mean(sigma_IVW_meta)
mean(sigma_IVWs)
mean(sigma_y_TwoStage)
mean(sigma_y_IVW)
mean(sigma_y_IVW_meta)

# n <- 15000
# p_thres = c(5E-04,1E-03,5E-03,1E-02,5E-02,1E-01,5E-01)
# n_thres <- length(p_thres)
# TwoStage_est = rep(0,times)
# IVW_est = rep(0,times)
# IVWs_est = rep(0,times)
# IVW_est1 = rep(0,times)
# cover_TwoStage_est = rep(0,times)
# cover_IVW_est = rep(0,times)
# cover_IVWs_est = rep(0,times)
# cover_IVW_est1 = rep(0,times)
# sigma_TwoStage = rep(0,times)
# sigma_IVW = rep(0,times)
# sigma_IVWs = rep(0,times)
# sigma_IVW1 = rep(0,times)
# 
# TwoStage_est_all = matrix(0,times,n_thres)
# IVW_est_all = matrix(0,times,n_thres)
# IVWs_est_all = matrix(0,times,n_thres)
# IVW_est_all1 = matrix(0,times,n_thres)


#sigma_est = rep(0,times)
#sigma_beta_est = rep(0,times)
# n.test = 1000
# n.snp = 100
# G_ori = matrix(rbinom((n+n.test)*n.snp,1,MAF),n+n.test,n.snp)
# G_all = apply(G_ori,2,scale)
# twostage.nsnps <- rep(0,times)
# twostage.prop <- rep(0,times)
# IVW.nsnps <- rep(0,times)
# IVW.prop <- rep(0,times)
# for(i in 1:times){
#   print(i)
#   beta_M = 0.1
#   beta_U = 0.1
#   alpha_G = 0.01
#   alpha_U = 0.01
#   U = rnorm(n+n.test)
#   p = 5
#   MAF=0.25
#   beta_G = rep(0.01,p)
#   sigma_y = 1
#   sigma_x = 1
#   #M = G%*%beta_G+rnorm(n,sd = sqrt(sigma_y))
#   #Y = M*beta_M +rnorm(n,sd = sqrt(sigma_x))
#   
#   M_all = G_all[,1:p]%*%beta_G+U*alpha_U+rnorm(n+n.test,sd = sqrt(sigma_y))
#   Y_all = M_all*beta_M + U*beta_U+rnorm(n+n.test,sd = sqrt(sigma_x))
#   M = M_all[1:n,]
#   Y = Y_all[1:n]
#   G = G_all[1:n,]
#   Y_test = Y_all[(n+1):(n+n.test)]
#   G_test = G_all[(n+1):(n+n.test),]
#   
#   p_vec = rep(0,n.snp)
#   beta_vec = rep(0,n.snp)
#   for(j in 1:n.snp){
#     model_temp = lm(M~G[,j]-1)
#     coef_temp = coef(summary(model_temp))
#     p_vec[j] <- coef_temp[4]
#   }
#  
#   r2_vec = rep(0,length(p_thres))
#   for(l in 1:length(p_thres)){
#     idx <- which(p_vec<=p_thres[l])
#     if(length(idx)==0){
#       r2_vec[l] = 0
#     }else{
#       model_temp2 = lm(Y~G[,idx]-1)
#       beta_temp = coefficients(model_temp2)
#       PRS_test = G_test[,idx,drop=F]%*%beta_temp
#       model_test = lm(Y_test~PRS_test-1)
#       r2_vec[l] = summary(model_test)$r.squared  
#     }
#     
#   }
#   
#   kdx = which(p_vec<=p_thres[which.max(r2_vec)])
#   twostage.nsnps[i] <- length(kdx)
#   twostage.prop[i] <- sum(c(1:5)%in%kdx)/p
# 
#     G_two = G[,kdx,drop=F]
#     
#     TwoStage_result = TwoStage(Y,M,G_two,beta_M)  
#     TwoStage_est[i] = TwoStage_result[1]
#     cover_TwoStage_est[i] = TwoStage_result[2]
#     sigma_TwoStage[i] = TwoStage_result[3]
#     
#   
#   
#   ldx = kdx
#   IVW.nsnps[i] <- length(ldx)
#   IVW.prop[i] <- sum(c(1:5)%in%ldx)/p
#   
#   if(length(ldx)==0){
#     IVW_result = NA
#     IVW_est[i] = NA
#     cover_IVW_est[i] = NA
#     sigma_IVW[i] = NA
#     IVWs_result = NA
#     IVWs_est[i] = NA
#     cover_IVWs_est[i] = NA
#     sigma_IVWs[i] = NA
#   }else{
#     G_oth = G[,ldx,drop=F]
#     IVW_result = IVW(Y,M,G_oth,beta_M)  
#     IVW_est[i] = IVW_result[1]
#     cover_IVW_est[i] = IVW_result[2]
#     sigma_IVW[i] = IVW_result[3]
#     IVW_est1[i] = IVW_result[5]
#     cover_IVW_est1[i] = IVW_result[6]
#     sigma_IVW1[i] = IVW_result[7]
#     
#     IVWs_result = IVW_s(Y,M,G_oth,beta_M)  
#     IVWs_est[i] = IVWs_result[1]
#     cover_IVWs_est[i] = IVWs_result[2]
#     sigma_IVWs[i] = IVWs_result[3]
#   }
#   
#   #get the estimate for every p.value threshold
#   for(l in 1:n_thres){
#     kdx = which(p_vec<=p_thres[l])
#     if(length(kdx)==0){
#       TwoStage_est_all[i,l] = NA
#       IVW_est_all[i,l] = NA
#       IVWs_est_all[i,l] = NA
#       IVW_est_all1[i,l] = NA
#     }else{
#       G_two = G[,kdx,drop=F]  
#       TwoStage_est_all[i,l] = TwoStage(Y,M,G_two,beta_M)[1]
#       IVW_est_all[i,l] = IVW(Y,M,G_two,beta_M)[1]
#       IVWs_est_all[i,l] = IVW_s(Y,M,G_oth,beta_M)[1]  
#       IVW_est_all1[i,l] = IVW(Y,M,G_two,beta_M)[5]
#     #twostage.nsnps[i] <- length(kdx)
#     #twostage.prop[i] <- sum(c(1:5)%in%kdx)/p
#      }
#   
#   
#   
#   
#   }
# }
# 
# 
# 
# result2 = list(TwoStage_est,IVW_est,IVWs_est,
#                IVW_est1,
#                cover_TwoStage_est,cover_IVW_est,
#                cover_IVWs_est,
#                cover_IVW_est1,
#                sigma_TwoStage,
#                sigma_IVW,
#                sigma_IVWs,
#                sigma_IVW1,
#                twostage.nsnps,
#                twostage.prop,
#                TwoStage_est_all,
#                IVW_est_all,
#                IVWs_est_all,
#                IVW_est_all1)


times = 100
n <- 150000
MAF =0.25
p <- 5
TwoStage_est = rep(0,times)
cover_TwoStage_est = rep(0,times)
sigma_TwoStage = rep(0,times)
sigma_y_TwoStage = rep(0,times)
IVW_est = rep(0,times)
cover_IVW_est = rep(0,times)
sigma_IVW = rep(0,times)
sigma_y_IVW = rep(0,times)
IVW_meta_est = rep(0,times)
cover_IVW_meta_est = rep(0,times)
sigma_IVW_meta = rep(0,times)
sigma_y_IVW_meta = rep(0,times)
IVWs_est = rep(0,times)
cover_IVWs_est = rep(0,times)
sigma_IVWs = rep(0,times)

#sigma_est = rep(0,times)
#sigma_beta_est = rep(0,times)
G_ori = matrix(rbinom(n*5,1,MAF),n,p)
G = apply(G_ori,2,scale)
for(i in 1:times){
  print(i)
  beta_M = 0.1
  beta_U = 0.1
  
  alpha_U = 0.01
  U = rnorm(n)
  # p = 5
  # MAF=0.25
  beta_G = rep(beta_G_vec[i2],p)
  sigma_y = 1
  sigma_m = 1
  M = G%*%beta_G+rnorm(n,sd = sqrt(sigma_m))
  Y = M*beta_M +rnorm(n,sd = sqrt(sigma_y))
  
  #M = G%*%beta_G+U*alpha_U+rnorm(n,sd = sqrt(sigma_m))
  #Y = M*beta_M + U*beta_U+rnorm(n,sd = sqrt(sigma_y))
  #M = G%*%beta_G+rnorm(n,sd = sqrt(sigma_m))
  #Y = M*beta_M + rnorm(n,sd = sqrt(sigma_y))
  TwoStage_result = TwoStage(Y,M,G,beta_M)  
  TwoStage_est[i] = TwoStage_result[1]
  cover_TwoStage_est[i] = TwoStage_result[2]
  sigma_TwoStage[i] = TwoStage_result[3]
  sigma_y_TwoStage[i] = TwoStage_result[4]
  IVW_result = IVW(Y,M,G,beta_M)  
  IVW_est[i] = IVW_result[1]
  cover_IVW_est[i] = IVW_result[2]
  sigma_IVW[i] = IVW_result[3]
  sigma_y_IVW[i] = IVW_result[4]
  IVW_meta_result = IVW_meta(Y,M,G,beta_M)  
  IVW_meta_est[i] = IVW_meta_result[1]
  cover_IVW_meta_est[i] = IVW_meta_result[2]
  sigma_IVW_meta[i] = IVW_meta_result[3]
  sigma_y_IVW_meta[i] = IVW_meta_result[4]
  IVWs_result = IVW_s(Y,M,G,beta_M)  
  IVWs_est[i] = IVWs_result[1]
  cover_IVWs_est[i] = IVWs_result[2]
  sigma_IVWs[i] = IVWs_result[3]
  
}

result3 = list(TwoStage_est,
               cover_TwoStage_est,
               sigma_TwoStage,
               sigma_y_TwoStage,
               IVW_est,
               cover_IVW_est,
               sigma_IVW,
               sigma_y_IVW,
               IVW_meta_est,
               cover_IVW_meta_est,
               sigma_IVW_meta,
               sigma_y_IVW_meta,
               IVWs_est,
               cover_IVWs_est,
               sigma_IVWs
)

# mean(TwoStage_est)-beta_M
# mean(IVW_est)-beta_M
# mean(IVWs_est)-beta_M
# mean(cover_TwoStage_est)
# mean(cover_IVW_est)
# mean(cover_IVWs_est,na.rm = T)
# var(TwoStage_est)
# var(IVW_est)
# var(IVWs_est)
# mean(sigma_TwoStage)
# mean(sigma_IVW)
# mean(sigma_IVWs)




# n <- 150000
# p_thres = c(5E-04,1E-03,5E-03,1E-02,5E-02,1E-01,5E-01)
# n_thres <- length(p_thres)
# TwoStage_est = rep(0,times)
# IVW_est = rep(0,times)
# IVWs_est = rep(0,times)
# IVW_est1 = rep(0,times)
# cover_TwoStage_est = rep(0,times)
# cover_IVW_est = rep(0,times)
# cover_IVWs_est = rep(0,times)
# cover_IVW_est1 = rep(0,times)
# sigma_TwoStage = rep(0,times)
# sigma_IVW = rep(0,times)
# sigma_IVWs = rep(0,times)
# sigma_IVW1 = rep(0,times)
# 
# TwoStage_est_all = matrix(0,times,n_thres)
# IVW_est_all = matrix(0,times,n_thres)
# IVWs_est_all = matrix(0,times,n_thres)
# IVW_est_all1 = matrix(0,times,n_thres)
# 
# 
# #sigma_est = rep(0,times)
# #sigma_beta_est = rep(0,times)
# n.test = 1000
# n.snp = 100
# G_ori = matrix(rbinom((n+n.test)*n.snp,1,MAF),n+n.test,n.snp)
# G_all = apply(G_ori,2,scale)
# twostage.nsnps <- rep(0,times)
# twostage.prop <- rep(0,times)
# IVW.nsnps <- rep(0,times)
# IVW.prop <- rep(0,times)
# for(i in 1:times){
#   print(i)
#   beta_M = 0.1
#   beta_U = 0.1
#   alpha_G = 0.01
#   alpha_U = 0.01
#   U = rnorm(n+n.test)
#   p = 5
#   MAF=0.25
#   beta_G = rep(0.01,p)
#   sigma_y = 1
#   sigma_x = 1
#   #M = G%*%beta_G+rnorm(n,sd = sqrt(sigma_y))
#   #Y = M*beta_M +rnorm(n,sd = sqrt(sigma_x))
#   
#   M_all = G_all[,1:p]%*%beta_G+U*alpha_U+rnorm(n+n.test,sd = sqrt(sigma_y))
#   Y_all = M_all*beta_M + U*beta_U+rnorm(n+n.test,sd = sqrt(sigma_x))
#   M = M_all[1:n,]
#   Y = Y_all[1:n]
#   G = G_all[1:n,]
#   Y_test = Y_all[(n+1):(n+n.test)]
#   G_test = G_all[(n+1):(n+n.test),]
#   
#   p_vec = rep(0,n.snp)
#   beta_vec = rep(0,n.snp)
#   for(j in 1:n.snp){
#     model_temp = lm(M~G[,j]-1)
#     coef_temp = coef(summary(model_temp))
#     p_vec[j] <- coef_temp[4]
#   }
#   
#   r2_vec = rep(0,length(p_thres))
#   for(l in 1:length(p_thres)){
#     idx <- which(p_vec<=p_thres[l])
#     if(length(idx)==0){
#       r2_vec[l] = 0
#     }else{
#       model_temp2 = lm(Y~G[,idx]-1)
#       beta_temp = coefficients(model_temp2)
#       PRS_test = G_test[,idx,drop=F]%*%beta_temp
#       model_test = lm(Y_test~PRS_test-1)
#       r2_vec[l] = summary(model_test)$r.squared  
#     }
#     
#   }
#   
#   kdx = which(p_vec<=p_thres[which.max(r2_vec)])
#   twostage.nsnps[i] <- length(kdx)
#   twostage.prop[i] <- sum(c(1:5)%in%kdx)/p
#   
#   G_two = G[,kdx,drop=F]
#   
#   TwoStage_result = TwoStage(Y,M,G_two,beta_M)  
#   TwoStage_est[i] = TwoStage_result[1]
#   cover_TwoStage_est[i] = TwoStage_result[2]
#   sigma_TwoStage[i] = TwoStage_result[3]
#   
#   
#   
#   ldx = kdx
#   IVW.nsnps[i] <- length(ldx)
#   IVW.prop[i] <- sum(c(1:5)%in%ldx)/p
#   
#   if(length(ldx)==0){
#     IVW_result = NA
#     IVW_est[i] = NA
#     cover_IVW_est[i] = NA
#     sigma_IVW[i] = NA
#     IVWs_result = NA
#     IVWs_est[i] = NA
#     cover_IVWs_est[i] = NA
#     sigma_IVWs[i] = NA
#   }else{
#     G_oth = G[,ldx,drop=F]
#     IVW_result = IVW(Y,M,G_oth,beta_M)  
#     IVW_est[i] = IVW_result[1]
#     cover_IVW_est[i] = IVW_result[2]
#     sigma_IVW[i] = IVW_result[3]
#     IVW_est1[i] = IVW_result[5]
#     cover_IVW_est1[i] = IVW_result[6]
#     sigma_IVW1[i] = IVW_result[7]
#     
#     IVWs_result = IVW_s(Y,M,G_oth,beta_M)  
#     IVWs_est[i] = IVWs_result[1]
#     cover_IVWs_est[i] = IVWs_result[2]
#     sigma_IVWs[i] = IVWs_result[3]
#   }
#   
#   #get the estimate for every p.value threshold
#   for(l in 1:n_thres){
#     kdx = which(p_vec<=p_thres[l])
#     if(length(kdx)==0){
#       TwoStage_est_all[i,l] = NA
#       IVW_est_all[i,l] = NA
#       IVWs_est_all[i,l] = NA
#       IVW_est_all1[i,l] = NA
#     }else{
#       G_two = G[,kdx,drop=F]  
#       TwoStage_est_all[i,l] = TwoStage(Y,M,G_two,beta_M)[1]
#       IVW_est_all[i,l] = IVW(Y,M,G_two,beta_M)[1]
#       IVWs_est_all[i,l] = IVW_s(Y,M,G_oth,beta_M)[1]  
#       IVW_est_all1[i,l] = IVW(Y,M,G_two,beta_M)[5]
#       #twostage.nsnps[i] <- length(kdx)
#       #twostage.prop[i] <- sum(c(1:5)%in%kdx)/p
#     }
#     
#     
#     
#     
#   }
# }
# 
# 
# 
# result4 = list(TwoStage_est,IVW_est,IVWs_est,
#                IVW_est1,
#                cover_TwoStage_est,cover_IVW_est,
#                cover_IVWs_est,
#                cover_IVW_est1,
#                sigma_TwoStage,
#                sigma_IVW,
#                sigma_IVWs,
#                sigma_IVW1,
#                twostage.nsnps,
#                twostage.prop,
#                TwoStage_est_all,
#                IVW_est_all,
#                IVWs_est_all,
#                IVW_est_all1)


result = list(result1,result3)
save(result,file = paste0("./result/simulation/simulation_M_fixed",i1,"_",i2,".Rdata"))



