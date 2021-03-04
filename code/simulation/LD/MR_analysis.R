args = commandArgs(trailingOnly = T)
l = as.numeric(args[[1]])
library(withr)
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/", install_github('qingyuanzhao/mr.raps'))
#install_github('qingyuanzhao/mr.raps')
library(mr.raps)
library(MendelianRandomization)
library(MRPRESSO)
setwd("/data/zhangh24/MR_MA/")
source("./code/simulation/functions/MR_function.R")
sum.data.y = as.data.frame(fread(paste0(cur.dir,"y_summary_chr_",j,"_rho_",l)))
sum.data.m = as.data.frame(fread(paste0(cur.dir,"m_summary_chr_",j,"_rho_",l)))
n.snp = nrow(sum.data.m)
n.rep = 100
for(i_rep in  1:n.rep){
  p = sum.data.m[,(6+3*i_rep)]
  print(range(p))
}
  idx = which(p<=0.05/n.snp)
  Gamma = sum.data.y[idx,(6+3*i_rep-2)]
  var_Gamma = (sum.data.y[idx,(6+3*i_rep-2)]/sum.data.y[idx,(6+3*i_rep-1)])^2
  gamma = sum.data.m[idx,(6+3*i_rep-2)]
  var_gamma = (sum.data.m[idx,(6+3*i_rep-2)]/sum.data.m[idx,(6+3*i_rep-1)])^2
  
  num = 3
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                      gamma,var_gamma)
  IVW_c_result[k] <- paste0(round(IVW_c_temp[[1]],num)," (",round(IVW_c_temp[[4]],num),")")
  MR_weight_temp = MRWeight(Gamma,var_Gamma,
                            gamma,var_gamma)
  keep.id = MR_weight_temp[[5]]
  out.id = MR_weight_temp[[6]]
  
  keep.ind = rep(FALSE,length(Gamma))
  keep.ind[out.id] = TRUE
  keep.ind = factor(keep.ind,levels=c(TRUE,FALSE))
  lm(Gamma[keep.id]~gamma[keep.id])
  
  beta_est = MR_weight_temp[[1]]
  plot.data = data.frame(gamma,Gamma,Weight = 1/(var_Gamma+beta_est^2*var_gamma),Removed = keep.ind)
  p = ggplot(plot.data,aes(gamma,Gamma))+
    geom_point(aes(col=Removed,size = Weight),alpha = 0.9)+
    geom_abline(slope = MR_weight_temp[[1]],intercept = 0)+
    theme_Publication()+
    scale_colour_Publication()+
    ggtitle(paste0(trait_vec[k]))+
    xlab("Alpha")+
    geom_abline(slope =-0.502,intercept = 0.006, linetype = "dashed")
  
  png(filename = paste0("./result/real_data_analysis/",trait_vec[k],"annotate.png"),width = 8,height = 6, units = "in",res= 300)
  print(p)
  dev.off()
  
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
                                           se.outcome = sqrt(var_Gamma)),
  )
  MR_raps_result[k] <- paste0(round(raps_result$beta.hat,num)," (",round(raps_result$beta.se,num),")")
  summary.data = data.frame(E1_effect = gamma,
                            E1_se = sqrt(var_gamma),
                            
                            Y_effect = Gamma,
                            Y_se = sqrt(var_Gamma))
  presso_result <- mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = summary.data, NbDistribution = 1000,  SignifThreshold = 0.05)
  mr_presso_result[k] = paste0(round(presso_result$`Main MR results`[1,3],num)," (",round(presso_result$`Main MR results`[1,4],num),")")
}