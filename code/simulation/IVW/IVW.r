#test the type one error and coverage prob of two-stage, IVW and IVW summary level statistics
#IVW estimate the variance using the meta-analysis variance
args = commandArgs(trailingOnly = T)
array = as.numeric(args[[1]])

array.mat <- matrix(0,3*4*4*100,4)
array.value <- rep(0,3*4*4*100)
temp <- 1
for(i1 in 1:3){
  for(i2 in 1:4){
    for(i3 in 1:100){
      for(i4 in 1:4){
        array.value[temp] <- 4*4*100*(i1-1)+4*100*(i2-1)+4*(i3-1)+i4
        array.mat[temp,1] = i1
        array.mat[temp,2] = i2
        array.mat[temp,3] = i3
        array.mat[temp,4] = i4
        temp = temp+1
      }
    }
  }
}
idx <- which(array.value==array)
i1 = array.mat[idx,1]
i2 = array.mat[idx,2]
i3 = array.mat[idx,3]
i4 = array.mat[idx,4]
# i1 = as.numeric(args[[1]])
# i2 = as.numeric(args[[2]])
# i3 = as.numeric(args[[3]])
# i4 = as.numeric(args[[4]])
print(i1)
print(i2)
print(i3)
print(i4)
setwd("/n/holystore01/LABS/xlin/Lab/hzhang/MR_MA")
#setwd("/data/zhangh24/MR_MA/")
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
  K <- length(Gamma)
  #return the plug in test results
  return(sum((Gamma-beta_plug*gamma)^2/(var_Gamma+beta_plug^2*var_gamma))-qchisq(0.95,K))
}

QuacNew <- function(Gamma,var_Gamma,gamma,var_gamma,beta_plug,W_vec){
  K <- length(Gamma)
  #return the plug in test results
  return(sum((Gamma-beta_plug*gamma)^2/(W_vec))-qchisq(0.95,K))
}

Gamma <- Gamma_est[idx[1],]
var_Gamma <- Gamma_var[idx[1],]
gamma <- Gamma_var[idx[1],]
var_gamma<- gamma_var[idx[1],]


ARMethod <- function(Gamma,var_Gamma,gamma,var_gamma){
  K <- length(Gamma)
  keep.ind <- c(1:K)

  beta_seq <- seq(-5,5,by=0.001)

  quan_result_AR <- rep(0,length(beta_seq))
  for(k in 1:length(beta_seq)){
     quan_result_AR[k] <- QuacForm(Gamma,var_Gamma,gamma,var_gamma,beta_seq[k])
  }
  #get the best estimate
  coef_best <- beta_seq[which.min(quan_result_AR)]
idx <- which(quan_result_AR<=0)
  if(length(idx)>0){
    beta_ci_range <- beta_seq[idx]
    jdx <- which(round(diff(beta_ci_range),3)!=0.001)
    if(length(jdx)==0){
      #remove the cases when AR potential wide confidence intervals
      coef_low <- min(beta_ci_range)
      coef_high <- max(beta_ci_range)
      cover_AR = ifelse((beta_M>=coef_low&
                           beta_M<=coef_high),1,0)
    }else{
      coef_low <- NA
      coef_high <- NA
      cover_AR = NA
      }
    
}else{
    coef_low = NA
    coef_high = NA
    cover_AR = 0
    
  }
  quan_result_true <- QuacForm(Gamma,var_Gamma,gamma,var_gamma,beta_M)
  cover_AR = ifelse(quan_result_true<=0,1,0)
  return(list(coef_best,coef_low,coef_high,cover_AR))
}


MRLR <- function(Gamma,var_Gamma,gamma,var_gamma){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  
  #first step
  model1 = lm(Gamma~gamma-1)
  coef_est = coefficients(model1)
  W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
 
  coef_best = sum(Gamma*gamma*W_vec)/sum(gamma^2*W_vec)
  sigma_est  = sum((Gamma-coef_est*gamma)^2)/(K-1)
  
  W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  xwx_iv = 1/sum(gamma^2*W_vec)
 
  var_coef_est = sigma_est*xwx_iv*t(gamma)%*%diag(W_vec)%*%diag(W_vec)%*%gamma*xwx_iv
 
  coef_low <- coef_est+qt(0.025,(K-1))*sqrt(var_coef_est)
  coef_high <- coef_est+qt(0.975,(K-1))*sqrt(var_coef_est)
  #coef_low_update <- confint(model1,level=0.95)[1]
  #coef_high_update <- confint(model1,level=0.95)[2]
  cover <- ifelse((beta_M>=coef_low&
                            beta_M<=coef_high),1,0)
 
  return(list(coef_est,coef_low,coef_high,cover))
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
ratio_AR_low <- rep(0,times)
ratio_AR_high <- rep(0,times)
ratio_est_MR <- rep(0,times)
ratio_MR_low <- rep(0,times)
ratio_MR_high <- rep(0,times)
cover_MR <- rep(0,times)

ind <- rep(0,times)
G_ori = matrix(rbinom(n*p,2,MAF),n,p)
G = apply(G_ori,2,scale)
G_ori2 = matrix(rbinom(n*p,2,MAF),n,p)
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
  ratio_est_c[i] <- as.numeric(ratio_temp[1])
  ratio_var_c[i] <- as.numeric(ratio_temp[2])
  ratio_cover_c[i] <- as.numeric(ratio_temp[3])
  ci_low_ratio_c[i] <- as.numeric(ratio_temp[4])
  ci_high_ratio_c[i] <- as.numeric(ratio_temp[5])
  
  ratio_AR_temp <- 
    ARMethod(as.numeric(est[[1]]),
             as.numeric(est[[2]]),
             as.numeric(est[[3]]),
             as.numeric(est[[4]]))
  ratio_est_AR[i] <- as.numeric(ratio_AR_temp[1])
  ratio_AR_low[i] <- as.numeric(ratio_AR_temp[2])
  ratio_AR_high[i] <- as.numeric(ratio_AR_temp[3])
  cover_AR[i] <- as.numeric(ratio_AR_temp[4])
  ratio_MR_temp <- 
    MRLR(as.numeric(est[[1]]),
             as.numeric(est[[2]]),
             as.numeric(est[[3]]),
             as.numeric(est[[4]]))
  ratio_est_MR[i] <- as.numeric(ratio_MR_temp[1])
  ratio_MR_low[i] <- as.numeric(ratio_MR_temp[2])
  ratio_MR_high[i] <- as.numeric(ratio_MR_temp[3])
  cover_MR[i] <- as.numeric(ratio_MR_temp[4])
  
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
              ratio_AR_low,
              ratio_AR_high,
              cover_AR,
              ratio_est_MR,
              ratio_MR_low,
              ratio_MR_high,
              cover_MR
)


save(result,file = paste0("./result/simulation/IVW/ratio_estimate_",i1,"_",i2,"_",i3,"_",i4,".Rdata"))



