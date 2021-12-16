#generate data based on summary statistics
n.snp = 13
h2 = 0.0367
N = 60000
n.rep = 100
beta_est_Raps = rep(NA,n.rep)
beta_cover_Raps = rep(NA,n.rep)
beta_se_Raps = rep(NA,n.rep)
beta_est_IVW = rep(NA,n.rep)
beta_cover_IVW = rep(NA,n.rep)
beta_se_IVW = rep(NA,n.rep)
beta_est_egger = rep(NA,n.rep)
beta_cover_egger = rep(NA,n.rep)
beta_se_egger  = rep(NA,n.rep)
beta_est_median = rep(NA,n.rep)
beta_cover_median = rep(NA,n.rep)
beta_se_median = rep(NA,n.rep)

library(withr)
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/", install_github('qingyuanzhao/mr.raps'))
#install_github('qingyuanzhao/mr.raps')
library(mr.raps)
library(MendelianRandomization)
library(MRPRESSO)
library(data.table)
library(dplyr)
library(MESS)
alpha = c(-0.158,
          0.129,
          0.271,-0.086,-0.2587,0.113,0.101,
          -0.055,-0.050,-0.044,0.0683,0.0345,-0.0599)
MAF_select = c(0.258,0.255,0.040,0.342,0.020,0.093,0.102,0.273,0.280,0.280,0.069,0.456,0.089)
for(i_rep in 1:n.rep){
  #generate alpha

  
  beta = 0
  scale.factor = sqrt(2*MAF_select*(1-MAF_select))
  alpha_select = rnorm(n.snp,alpha,sd = sqrt(1/N)/scale.factor)
  se_alpha_select = rep(sqrt(1/N),n.snp)
  se_alpha_select=  se_alpha_select/scale.factor
  
  
  Gamma_select = rnorm(n.snp,beta*alpha,sd = sqrt(1/N)/scale.factor)
  se_Gamma_select = rep(sqrt(1/N),n.snp)
 
  se_Gamma_select = se_Gamma_select/scale.factor
  MRInputObject <- mr_input(bx = alpha_select,
                            bxse = se_alpha_select,
                            by = Gamma_select,
                            byse = se_Gamma_select)
  
  IVWObject <- mr_ivw(MRInputObject,
                      model = "default",
                      robust = FALSE,
                      penalized = FALSE,
                      correl = FALSE,
                      weights = "simple",
                      psi = 0,
                      distribution =
                        "normal",
                      alpha = 0.05)
  
  
  beta_est_IVW[i_rep] = IVWObject$Estimate
  beta_cover_IVW[i_rep] = ifelse(IVWObject$CILower<=beta&
                                   IVWObject$CIUpper>=beta,1,0)
  beta_se_IVW[i_rep] = IVWObject$StdError
  EggerObject <- mr_egger(
    MRInputObject,
    robust = FALSE,
    penalized = FALSE,
    correl = FALSE,
    distribution = "normal",
    alpha = 0.05
  )
  
  beta_est_egger[i_rep] = EggerObject$Estimate
  beta_cover_egger[i_rep] = ifelse(EggerObject$CILower.Est<=beta&
                                     EggerObject$CIUpper.Est>=beta,1,0)
  beta_se_egger[i_rep] = EggerObject$StdError.Est
  
  MedianObject <- mr_median(
    MRInputObject,
    weighting = "weighted",
    distribution = "normal",
    alpha = 0.05,
    iterations = 10000,
    seed = 314159265
  )
  
  
  beta_est_median[i_rep] = MedianObject$Estimate
  beta_cover_median[i_rep] = ifelse(MedianObject$CILower<=beta&
                                      MedianObject$CIUpper>=beta,1,0)
  beta_se_median[i_rep] = MedianObject$StdError
  
  
  raps_result <- mr.raps(data = data.frame(beta.exposure = alpha_select,
                                           beta.outcome = Gamma_select,
                                           se.exposure = se_alpha_select,
                                           se.outcome = se_Gamma_select),
                         diagnostics = F)
  beta_est_Raps[i_rep] = raps_result$beta.hat
  beta_cover_Raps[i_rep] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                                    raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  beta_se_Raps[i_rep] = raps_result$beta.se  
}
method = c("IVW","MR-Egger","MR-median","MRRAPs")

mean.result = data.frame(
  beta_est_IVW,beta_est_egger,beta_est_median,beta_est_Raps
)
colnames(mean.result) = method

se.result = data.frame(
  beta_se_IVW,beta_se_egger,beta_se_median,beta_se_Raps
)
colnames(se.result) = method

cover.result = data.frame(
  beta_cover_IVW,beta_cover_egger,beta_cover_median,beta_cover_Raps
)
colnames(cover.result) = method

result = data.frame(
  method = method,
  bias = apply(mean.result,2,function(x){mean(x,na.rm=T)})-beta,
  em_se = apply(mean.result,2,function(x){sd(x,na.rm=T)}),
  es_se = apply(se.result,2,function(x){mean(x,na.rm=T)}),
  cover = apply(cover.result,2,function(x){mean(x,na.rm=T)}),
  rmse = apply(mean.result,2,function(x){sqrt(mean((x-beta)^2,na.rm=T))})
)
print(result)
result$i_vec = rep(i,length(method))
result$l_vec = rep(l,length(method))
result$v_vec = rep(v,length(method))
