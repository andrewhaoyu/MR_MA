#bmi analysis 
#first run LD clumping to get genome-wide significant SNPs
#second run the IVW and beta estimate model
#update date: 012720
#down load data ( https://github.com/qingyuanzhao/mr.raps)

#setwd("/data/zhangh24/MR_MA")
library(data.table)
library(dplyr)
library(tidyr)
library(devtools)
library(withr)
library(mr.raps)
pcut <- c(1)
l <- length(pcut)
IVW_c_result <- rep("c",l)
MRweight_result <- rep("c",l)
MR_egger_result <- rep("c",l)
MR_median_result <- rep("c",l)
MR_raps_result <- rep("c",l)
MR_presso_result <- rep("c",l)
n.snp <- rep(0,l)
keep.snp <- rep(0,l)
num = 3
for(k in 1:l){
  pdx <- which(cad.cad$pval.selection<=pcut[k])
  out_clump_SNP_temp = cad.cad[pdx,]
  Gamma = out_clump_SNP_temp$beta.outcome
  var_Gamma = out_clump_SNP_temp$se.outcome^2
  gamma = out_clump_SNP_temp$beta.exposure
  var_gamma = out_clump_SNP_temp$se.exposure^2
  n.snp[k] <- length(Gamma)
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                      gamma,var_gamma)
  IVW_c_result[k] <- paste0(round(IVW_c_temp[[1]],num)," (",round(IVW_c_temp[[4]],num),")")
  MR_weight_temp = MRWeight(Gamma,var_Gamma,
                            gamma,var_gamma)
  MRweight_result[k] <- paste0(round(MR_weight_temp[[1]],num)," (",round(MR_weight_temp[[4]],num),")")
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
  MR_median_result[k] <- paste0(round(median_result$Estimate,num)," (",round(median_result$StdError,num),")")
  egger_result <- mr_egger(MRInputObject,
                           robust = FALSE,
                           penalized = FALSE,
                           correl = FALSE,
                           distribution = "normal",
                           alpha = 0.05)
  MR_egger_result[k] <- paste0(round(egger_result$Estimate,num)," (",round(egger_result$StdError.Est,num),")")
  raps_result <- mr.raps(data = data.frame(beta.exposure = gamma,
                                           beta.outcome = Gamma,
                                           se.exposure = sqrt(var_gamma),
                                           se.outcome = sqrt(var_Gamma)))
  MR_raps_result[k] <- paste0(round(raps_result$beta.hat,num)," (",round(raps_result$beta.se,num),")")
  summary.data = data.frame(E1_effect = gamma,
                            E1_se = sqrt(var_gamma),
                            
                            Y_effect = Gamma,
                            Y_se = sqrt(var_Gamma))
  presso_result <- mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = summary.data, NbDistribution = 1000,  SignifThreshold = 0.05)
  mr_presso_result = paste0(round(presso_result$`Main MR results`[1,3],num)," (",round(presso_result$`Main MR results`[1,4],num),")")
  # 
}
cad.result.summary <- data.frame(IVW_c_result,MRweight_result,MR_egger_result,MR_median_result
                                 ,MR_raps_result,mr_presso_result)
