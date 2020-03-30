#test the type one error and coverage prob of two-stage, IVW and IVW summary level statistics
#IVW estimate the variance using the meta-analysis variance
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])
i3 = as.numeric(args[[3]])
i4 = as.numeric(args[[4]])


setwd("/data/zhangh24/MR_MA/")
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





Regression = function(Y,M,G,G2){
  p <- ncol(G)
  Gamma <- rep(0,p)
  var_Gamma <- rep(0,p)
  gamma <- rep(0,p)
  var_gamma <- rep(0,p)
  for(k in 1:p){
    n = length(Y)
    model1 = lm(Y~G[,k]-1)
    coef_temp = coef(summary(model1))
    Gamma[k] = coef_temp[1]
    var_Gamma[k] = coef_temp[2]^2
    model2 = lm(M~G2[,k]-1)
    coef_temp2 = coef(summary(model2))
    gamma[k] = coef_temp2[1]
    var_gamma[k] = coef_temp2[2]^2
    
  }
  return(list(Gamma,var_Gamma,
           gamma,var_gamma))
  
  
}
#IVW estimate using summary level statistics
IVW_s = function(Gamma,var_Gamma,gamma,var_gamma){
  p <- length(Gamma)
  raio_vec = rep(0,p)
  ratio_var_vec = rep(0,p)
  
    raio_vec = Gamma/gamma
    ratio_var_vec =  var_Gamma/gamma^2
  
  
  Meta_result = Meta(raio_vec,ratio_var_vec)
  ratio_ivw =   Meta_result[1]
  ratio_ivw_var = Meta_result[2]
  coef_low = ratio_ivw-1.96*sqrt(ratio_ivw_var)
  coef_high = ratio_ivw+1.96*sqrt(ratio_ivw_var)
  cover = ifelse((beta_M>=coef_low&
                    beta_M<=coef_high),1,0)
  
  return(c(ratio_ivw,ratio_ivw_var,cover,
           coef_low,coef_high))
}

#IVW estimate using summary level statistics
IVW_c = function(Gamma,var_Gamma,gamma,var_gamma){
  p <- length(Gamma)
  raio_vec = rep(0,p)
  ratio_var_vec = rep(0,p)
  
  raio_vec = Gamma/gamma
  ratio_var_vec =  var_Gamma/gamma^2+var_gamma*Gamma^2/gamma^4
  
  
  Meta_result = Meta(raio_vec,ratio_var_vec)
  ratio_ivw =   Meta_result[1]
  ratio_ivw_var = Meta_result[2]
  coef_low = ratio_ivw-1.96*sqrt(ratio_ivw_var)
  coef_high = ratio_ivw+1.96*sqrt(ratio_ivw_var)
  cover = ifelse((beta_M>=coef_low&
                    beta_M<=coef_high),1,0)
  
  return(c(ratio_ivw,ratio_ivw_var,cover,
           coef_low,coef_high))
}



QuacForm <- function(Gamma,var_Gamma,gamma,var_gamma,beta_plug){
  return(sum((Gamma-beta_plug*gamma)^2/(var_Gamma+beta_plug^2*var_gamma)))
}







ARMethod <- function(Gamma,var_Gamma,gamma,var_gamma){
  beta_seq <- seq(-10,10,by=0.01)
  quan_result <- rep(0,length(beta_seq))
  for(k in 1:length(beta_seq)){
    quan_result[k] <- QuacForm(Gamma,var_Gamma,gamma,var_gamma,beta_seq[k])
  }
  coef_est <- beta_seq[which.min(quan_result)]
  test_result <- QuacForm(Gamma,var_Gamma,gamma,var_gamma,beta_M)
  
  cover = ifelse(test_result<=qchisq(0.95,df=length(Gamma)),1,0)
  
  return(c(coef_est,cover))
}

n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.00,0.01,0.03,0.05)
beta_vec = c(0,0.3,0.5,1)
times = 1000
p <- 5
n <- n_vec[i1]
MAF =0.25
Gamma_est <- matrix(0,times,p)
Gamma_var <- matrix(0,times,p)
gamma_est <- matrix(0,times,p)
gamma_var <- matrix(0,times,p)
ratio_est <- rep(0,times)
ratio_var <- rep(0,times)
ratio_cover <- rep(0,times)
ci_low_ratio<- rep(0,times)
ci_high_ratio <- rep(0,times)
ratio_est_c <- rep(0,times)
ratio_var_c <- rep(0,times)
ratio_cover_c <- rep(0,times)
ratio_cover_c <- rep(0,times)
ci_low_ratio_c <- rep(0,times)
ci_high_ratio_c <- rep(0,times)

ratio_est_AR <- rep(0,times)
cover_AR <- rep(0,times)
ind <- rep(0,times)
G_ori = matrix(rbinom(n*5,2,MAF),n,5)
G = apply(G_ori,2,scale)
G_ori2 = matrix(rbinom(n*5,2,MAF),n,5)
G2 = apply(G_ori2,2,scale)
set.seed(i3)
for(i in 1:times){
  print(i)
  beta_M = beta_vec[i4]
  alpha_G = rep(alpha_vec[i2],p)
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
  
  Gamma_est[i,] <- as.numeric(est[[1]])
  Gamma_var[i,] <- as.numeric(est[[2]])
  gamma_est[i,] <- as.numeric(est[[3]])
  gamma_var[i,] <- as.numeric(est[[4]])
  ratio_temp <- IVW_s(as.numeric(est[[1]]),
                      as.numeric(est[[2]]),
                      as.numeric(est[[3]]),
                      as.numeric(est[[4]]))
  ratio_est[i] <- ratio_temp[1]
  ratio_var[i] <- ratio_temp[2]
  ratio_cover[i] <- ratio_temp[3]
  ci_low_ratio[i] <- ratio_temp[4]
  ci_high_ratio[i] <- ratio_temp[5]
  ratio_temp <- IVW_c(as.numeric(est[[1]]),
                      as.numeric(est[[2]]),
                      as.numeric(est[[3]]),
                      as.numeric(est[[4]]))
  ratio_est_c[i] <- ratio_temp[1]
  ratio_var_c[i] <- ratio_temp[2]
  ratio_cover_c[i] <- ratio_temp[3]
  ci_low_ratio_c[i] <- ratio_temp[4]
  ci_high_ratio_c[i] <- ratio_temp[5]
  
  ratio_exact_temp <- 
    ARMethod(as.numeric(est[[1]]),
             as.numeric(est[[2]]),
             as.numeric(est[[3]]),
             as.numeric(est[[4]]))
  ratio_est_AR[i] <- ratio_exact_temp[1]
  cover_AR[i] <- ratio_exact_temp[2]
}




result = list(Gamma_est,
               Gamma_var,
               gamma_est,
               gamma_var,
               ratio_est,
               ratio_var,
               ratio_cover,
               ci_low_ratio,
               ci_high_ratio,
               ratio_est_c,
               ratio_var_c,
               ratio_cover_c,
               ratio_cover_c,
               ci_low_ratio_c,
               ci_high_ratio_c,
               ratio_est_AR,
               cover_AR
)


save(result,file = paste0("./result/simulation/IVW/ratio_estimate_",i1,"_",i2,"_",i3,"_",i4,".Rdata"))



