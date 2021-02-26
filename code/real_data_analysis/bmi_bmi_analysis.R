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
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/", install_github('qingyuanzhao/mr.raps'))
#install_github('qingyuanzhao/mr.raps')
library(mr.raps)
library(MendelianRandomization)
library(MRPRESSO)
#library(mr.raps, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/")
idx <- which(bmi.bmi$pval.selection<=5E-08)
Gamma = bmi.bmi$beta.outcome
var_Gamma = bmi.bmi$se.outcome^2
gamma = bmi.bmi$beta.exposure
var_gamma = bmi.bmi$se.exposure^2

pcut <- c(5E-08)
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
  pdx <- which(bmi.bmi$pval.selection<=pcut[k])
  out_clump_SNP_temp = bmi.bmi[pdx,]
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

bmi.result.summary <- data.frame(IVW_c_result,MRweight_result,MR_egger_result,MR_median_result
                                 ,MR_raps_result,mr_presso_result)
write.csv(bmi.result.summary,file = "/data/zhangh24/MR_MA/result/real_data_analysis/bmi/bmi.bmi.summary.csv")

lm(Gamma~gamma)


MRLR <- function(Gamma,var_Gamma,gamma,var_gamma){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  
  #first step
  model1 = lm(Gamma~gamma-1)
  coef_est = coefficients(model1)
  W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  
  coef_best = sum(Gamma*gamma*W_vec)/sum(gamma^2*W_vec)
  
  coef_vec <- seq(coef_best-0.1,coef_best+0.1,0.001)
  
  quad_vec = rep(0,length(coef_vec))
  for(i in 1:length(coef_vec)){
    coef_temp = coef_vec[i]
    W_vec = 1/(var_Gamma+coef_temp^2*var_gamma)
    quad_vec[i] = sum((Gamma-coef_temp*gamma)^2*W_vec)
  }
  
  coef_best = coef_vec[which.min(quad_vec)]
  return(coef_best)
}
MRLRVariance <- function(Gamma,var_Gamma,gamma,var_gamma,coef_est){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  
  sigma_est  = sum((Gamma-coef_est*gamma)^2)/(K-1)
  
  W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  xwx_iv = 1/sum(gamma^2*W_vec)
  var_coef_est = sigma_est*xwx_iv*t(gamma)%*%diag(W_vec)%*%diag(W_vec)%*%gamma*xwx_iv
  coef_low <- coef_est+qt(0.025,(K-1))*sqrt(var_coef_est)
  coef_high <- coef_est+qt(0.975,(K-1))*sqrt(var_coef_est)
  return(c(var_coef_est,
         coef_low,
         coef_high))
}


TestPleotropic <- function(Gamma,var_Gamma,gamma,var_gamma,coef_best){
  beta_plug = coef_best
  
  Gamma_update = Gamma
  var_Gamma_update = var_Gamma
  gamma_update = gamma
  var_gamma_update = var_gamma
  p_test_value = pchisq((Gamma_update-beta_plug*gamma_update)^2/(var_Gamma_update+beta_plug^2*var_gamma_update),1,lower.tail = F)
  p_adjust_test_value = p.adjust(p_test_value,method="none")  
  return(p_adjust_test_value)
}







MRWeight <- function(Gamma,var_Gamma,gamma,var_gamma){
  coef_best = MRLR(Gamma,var_Gamma,gamma,var_gamma)
  #test pleotropic
  K <- length(Gamma)
  all.id <- c(1:K)
  #id of SNPs that are removed
  out.id <- NULL
  #id of SNPs that are kept
  keep.id <- all.id
 
  p_adjust_test_value <- TestPleotropic(Gamma,var_Gamma,gamma,var_gamma,coef_best)
  cut.off = 0.01
  if(min(p_adjust_test_value)<cut.off){
  while(min(p_adjust_test_value)<cut.off){
    
    out.id <- c(out.id,keep.id[which.min(p_adjust_test_value)])
    keep.id <- setdiff(all.id,out.id)
    Gamma.keep = Gamma[keep.id]
    var_Gamma.keep = var_Gamma[keep.id]
    gamma.keep = gamma[keep.id]
    var_gamma.keep = var_gamma[keep.id]
    coef_best = MRLR(Gamma.keep,
                     var_Gamma.keep,
                     gamma.keep,
                     var_gamma.keep)
    p_adjust_test_value <- TestPleotropic(Gamma.keep,
                                          var_Gamma.keep,
                                          gamma.keep,
                                          var_gamma.keep,
                                          coef_best)
  }
  
  
}else{
  Gamma.keep = Gamma[keep.id]
  var_Gamma.keep = var_Gamma[keep.id]
  gamma.keep = gamma[keep.id]
  var_gamma.keep = var_gamma[keep.id]
}
  coef_est = coef_best
  result.temp <- MRLRVariance(Gamma.keep,
               var_Gamma.keep,
               gamma.keep,
               var_gamma.keep,
               coef_best)
  coef_var = result.temp[1]
  coef_low = result.temp[2]
  coef_high = result.temp[3]
  #coef_low_update <- confint(model1,level=0.95)[1]
  #coef_high_update <- confint(model1,level=0.95)[2]
  # cover <- ifelse((beta_M>=coef_low&
  #                    beta_M<=coef_high),1,0)
  # 
  return(list(coef_est,coef_low,coef_high,sqrt(coef_var),keep.id, out.id))
}


IVW_c <- function(Gamma,var_Gamma,gamma,var_gamma){
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
  
  
  return(c(ratio_ivw,
           coef_low,coef_high,sqrt(ratio_ivw_var)))
}
















sig_SNPs <- data.frame(out_clump_SNP$SNP,
                       as.numeric(out_clump_SNP$CHR),
                       as.numeric(out_clump_SNP$BP),
                       as.numeric(out_clump_SNP$P),
                       stringsAsFactors = F)
colnames(sig_SNPs) <- c("SNP","CHR","position",
                        "p.value")
PositionPruning <- function(sig_SNPs){
  sig_SNPs_temp =sig_SNPs
  filter_result = NULL
  temp.ind = 1
  while(nrow(sig_SNPs_temp)!=0){
    
    idx = which.min(sig_SNPs_temp$p.value)
    filter_result = rbind(filter_result,sig_SNPs_temp[idx,])
    
    position.range <- 500*10^3
    filter_result_position = sig_SNPs_temp$position[idx]
    filter_CHR = sig_SNPs_temp$CHR[idx]
    idx.cut <- which((sig_SNPs_temp$position>=filter_result_position-position.range)&(sig_SNPs_temp$position<=filter_result_position+position.range)&(sig_SNPs_temp$CHR==filter_CHR ))
    
    
    sig_SNPs_temp = sig_SNPs_temp[-idx.cut,]
    temp.ind = temp.ind+1
  }
  return(filter_result)
}

out_clump_SNP_clean <- PositionPruning(sig_SNPs)[,1,drop=F]
out_clump_SNP_clean <- left_join(out_clump_SNP_clean,
                                 out_clump_SNP,
                                 by="SNP")

#order the SNPs by chr and position
idx.order <- order(out_clump_SNP_clean$CHR,
                   out_clump_SNP_clean$BP)
out_clump_SNP_clean = out_clump_SNP_clean[idx.order,]
