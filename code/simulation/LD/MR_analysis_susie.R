args = commandArgs(trailingOnly = T)
l = as.numeric(args[[1]])
sub = as.numeric(args[[2]])
library(withr)
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/", install_github('qingyuanzhao/mr.raps'))
#install_github('qingyuanzhao/mr.raps')
library(mr.raps)
library(MendelianRandomization)
library(MRPRESSO)
library(data.table)
library(dplyr)
setwd("/data/zhangh24/MR_MA/")
source("./code/simulation/functions/MR_function.R")
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
j =22
sum.data.y = as.data.frame(fread(paste0(cur.dir,"y_summary_chr_",j,"_rho_",l)))
sum.data.m = as.data.frame(fread(paste0(cur.dir,"m_summary_chr_",j,"_rho_",l)))
pos = sum.data.m$BP
n.snp = nrow(sum.data.m)
n.rep = 100
library(susieR)
library(caret)
num = 10
library(bc2)
start.end = startend(n.rep,num,sub)
start = start.end[1]
end = start.end[2]
total = (end-start+1) *6
MR_result <- matrix(NA,total,4)

colnames(MR_result) <- c("est","sd","cover","method")
temp = 1
beta_M = 0.15
load(paste0(cur.dir,"chr_",j,"_LDmat.rdata"))
R = as.matrix(corr0)
for(i_rep in  start:end){
  Z = sum.data.m[,(6+3*i_rep-1)]
  
  p = sum.data.m[,(6+3*i_rep)]
  idx = which(p<=1E-03)
  for(k in 1:length(idx)){
    select.bp = pos[idx[k]]
    select.idx = which(pos>=select.bp-500000&pos<=select.bp+500000)
    z.select = Z[select.idx]
    R.select = R[select.idx,select.idx]
    #det(R.select[1:63,1:63])
    R.select.sub = R.select[-drop,-drop]
    z.select.sub = z.select[-drop] 
    fit_rss  = susie_rss(z.select.sub, R.select.sub, L = 10,
                         estimate_residual_variance = TRUE, 
                         estimate_prior_variance = TRUE)
    
  }
  
  
  
  fit_rss  = susie_rss(Z, R, L = 10,
            estimate_residual_variance = TRUE, 
            estimate_prior_variance = TRUE)
  idx = which(p<=0.05/n.snp)
  if(length(idx)>1){
    Gamma = sum.data.y[idx,(6+3*i_rep-2)]
    var_Gamma = (sum.data.y[idx,(6+3*i_rep-2)]/sum.data.y[idx,(6+3*i_rep-1)])^2
    gamma = sum.data.m[idx,(6+3*i_rep-2)]
    var_gamma = (sum.data.m[idx,(6+3*i_rep-2)]/sum.data.m[idx,(6+3*i_rep-1)])^2
    
    num = 3
    IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                        gamma,var_gamma)
    MR_result[temp,1] = c(IVW_c_temp[1])
    MR_result[temp,2] = c(IVW_c_temp[4])
    MR_result[temp,3] = ifelse(IVW_c_temp[2]<=beta_M&
                                 IVW_c_temp[3]>=beta_M,1,0)
    MR_result[temp,4] = "IVW"
    temp   = temp + 1
    MR_weight_temp = MRWeight(Gamma,var_Gamma,
                              gamma,var_gamma)
    MR_result[temp,1] = as.numeric(MR_weight_temp[[1]][1])
    MR_result[temp,2] = c(MR_weight_temp[[4]][1])
    MR_result[temp,3] = ifelse(MR_weight_temp[[2]][1]<=beta_M&
                                 MR_weight_temp[[3]][1]>=beta_M,1,0)
    MR_result[temp,4] = "MR-Weight"
    temp = temp+1
    MRInputObject <- mr_input(bx = gamma,
                              bxse = sqrt(var_gamma),
                              by = Gamma,
                              byse = sqrt(var_Gamma))
    median_result <- mr_median(MRInputObject,
                               weighting = "weighted",
                               distribution = "normal",
                               alpha = 0.05,
                               iterations = 10000,
                               seed = 314159265)
    MR_result[temp,1] = median_result$Estimate
    MR_result[temp,2] = median_result$StdError
    MR_result[temp,3] = ifelse(median_result$CILower<=beta_M&
                                 median_result$CIUpper>=beta_M,1,0)
    MR_result[temp,4] = "MR-Median"
    temp = temp +1
    egger_result <- mr_egger(MRInputObject,
                             robust = FALSE,
                             penalized = FALSE,
                             correl = FALSE,
                             distribution = "normal",
                             alpha = 0.05)
    MR_result[temp,1] = egger_result$Estimate
    MR_result[temp,2] = egger_result$StdError.Est
    MR_result[temp,3] = ifelse(egger_result$CILower.Est<=beta_M&
                                 egger_result$CIUpper.Est>=beta_M,1,0)
    MR_result[temp,4] = "MR-Egger"
    temp = temp +1
    raps_result <- mr.raps(data = data.frame(beta.exposure = gamma,
                                             beta.outcome = Gamma,
                                             se.exposure = sqrt(var_gamma),
                                             se.outcome = sqrt(var_Gamma)),
                           diagnostics = F)
    MR_result[temp,1] = raps_result$beta.hat
    MR_result[temp,2] = raps_result$beta.se
    MR_result[temp,3] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta_M&
                                 raps_result$beta.hat+1.96*raps_result$beta.se>=beta_M,1,0)
    MR_result[temp,4] = "Raps"
    temp = temp+1
    summary.data = data.frame(E1_effect = gamma,
                              E1_se = sqrt(var_gamma),
                              
                              Y_effect = Gamma,
                              Y_se = sqrt(var_Gamma))
    presso_result <- mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = summary.data, NbDistribution = 1000,  SignifThreshold = 0.05)
    MR_result[temp,1] = presso_est = presso_result$`Main MR results`[1,3]
    MR_result[temp,2] = presso_sd = presso_result$`Main MR results`[1,4]
    MR_result[temp,3] = ifelse(presso_est-1.96*presso_sd<=beta_M&
                                 presso_est+1.96*presso_sd>=beta_M,1,0)
    MR_result[temp,4] = "MR-PRESSO"
    
  }else{
    MR_result[temp,4] = "IVW"
    temp = temp+1
    MR_result[temp,4] = "MR-Weight" 
    temp  = temp+1
    MR_result[temp,4] = "MR-Median"
    temp = temp+1
    MR_result[temp,4] = "MR-Egger"
    temp = temp+1
    MR_result[temp,4] = "Raps"
    temp = temp+1
    MR_result[temp,4] = "MR-PRESSO"
    temp = temp+1
    
  }
}
save(MR_result,file =paste0(cur.dir,"MR_result_rho_",l,"_sub_",sub,".rdata"))
