#bmi analysis 
#first run LD clumping to get genome-wide significant SNPs
#second run the IVW and beta estimate model
#update date: 012720
#down load data ( https://github.com/qingyuanzhao/mr.raps)

setwd("/data/zhangh24/MR_MA")
library(data.table)
library(dplyr)
library(tidyr)
library(devtools)
library(withr)
with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/", install_github('qingyuanzhao/mr.raps'))
library(mr.raps, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
cad.cad <- cad.cad[-984,]
pcut <- c(5E-08,5E-07,5E-6,5E-5,5E-04,5E-03,5E-02,5E-01,1)
l <- length(pcut)
IVW_s_result <- rep("c",l)
IVW_c_result <- rep("c",l)
AR_result_1 <- rep("c",l)
AR_result_2 <- rep("c",l)
n.snp <- rep(0,l)
keep.snp <- rep(0,l)
for(k in 1:l){
  pdx <- which(cad.cad$pval.exposure<=pcut[k])
  out_clump_SNP_temp = cad.cad[pdx,]
  Gamma = out_clump_SNP_temp$beta.outcome
  var_Gamma = out_clump_SNP_temp$se.outcome^2
  gamma = out_clump_SNP_temp$beta.exposure
  var_gamma = out_clump_SNP_temp$se.exposure^2
  n.snp[k] <- length(Gamma)
  IVW_s_temp <- IVW_s(Gamma,var_Gamma,
                      gamma,var_gamma)
  num = 3
  IVW_s_result[k] <- paste0(round(IVW_s_temp[[1]],num)," (",round(IVW_s_temp[[3]],num),",",
                            round(IVW_s_temp[[4]],num),")")
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                      gamma,var_gamma)
  IVW_c_result[k] <- paste0(round(IVW_c_temp[[1]],num)," (",round(IVW_c_temp[[3]],num),",",
                            round(IVW_c_temp[[4]],num),")")
  AR_result <- ARMethod(Gamma,var_Gamma,
                        gamma,var_gamma)
  keep.snp[k] <- length(AR_result[[4]])
  AR_result_1[k] <- paste0(round(AR_result[[1]],num)," (",round(AR_result[[2]],num),",",
                           round(AR_result[[3]],num),")")
  AR_result_2[k] <- paste0(round(AR_result[[1]],num)," (",round(AR_result[[6]],num),",",
                           round(AR_result[[7]],num),")")
}
cad.result.summary <- data.frame(pcut,n.snp,IVW_s_result,
                                 IVW_c_result,AR_result_1,AR_result_2,keep.snp,stringsAsFactors = F)
write.csv(cad.result.summary,file = "/data/zhangh24/MR_MA/result/real_data_analysis/cad.cad.summary.csv")

