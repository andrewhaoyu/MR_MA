args = commandArgs(trailingOnly = T)
l = as.numeric(args[[1]])
v = as.numeric(args[[2]])
i = as.numeric(args[[3]])
#sub = as.numeric(args[[3]])
library(withr)
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/", install_github('qingyuanzhao/mr.raps'))
#install_github('qingyuanzhao/mr.raps')
library(mr.raps)
library(MendelianRandomization)
library(MRPRESSO)
library(data.table)
library(dplyr)
library(MESS)


temp = 1
result.list = list()
# for(l in 1:3){
#   for(v in 1:3){
   # for(i in 1:2){
      beta_vec = c(0,0.2)
      beta = beta_vec[i]
      #   for(sub in 1:10){
      setwd("/data/zhangh24/MR_MA/")
      source("./code/simulation/functions/MR_function.R")
      cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
      j =22
      sum.data.y = as.data.frame(fread(paste0(cur.dir,"y_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v)))
      sum.data.m = as.data.frame(fread(paste0(cur.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v)))
      n.snp = nrow(sum.data.m)
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
      for(i_rep in  1:n.rep){
        
        LD.snp = as.data.frame(fread( paste0(cur.dir,"LD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep,".clumped")))
        sum.data.match.m = left_join(LD.snp,sum.data.m,by=c("SNP"="ID"))
        p = sum.data.match.m[,(6+3*i_rep)]
        
        #idx = which(p<=0.05/nrow(sum.data.m))
        idx = which(p<=5E-08)
        if(length(idx)>3){
          sum.data.match.y = left_join(LD.snp,sum.data.y,by=c("SNP"="ID"))
          Gamma = sum.data.match.y[,(6+3*i_rep-2)]
          se_Gamma = as.numeric(sum.data.match.y[,(6+3*i_rep-1)])
          alpha = as.numeric(sum.data.match.m[,(6+3*i_rep-2)])
          se_alpha = as.numeric(sum.data.match.m[,(6+3*i_rep-1)])
          
          
          alpha_select =alpha[idx]
          se_alpha_select = se_alpha[idx]
          Gamma_select = Gamma[idx]
          se_Gamma_select = se_Gamma[idx]
          
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
          # summary.data = data.frame(E1_effect = gamma,
          #                           E1_se = sqrt(var_gamma),
          #                           
          #                           Y_effect = Gamma,
          #                           Y_se = sqrt(var_Gamma))
          # presso_result <- mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = summary.data, NbDistribution = 1000,  SignifThreshold = 0.05)
          # MR_result[temp,1] = presso_est = presso_result$`Main MR results`[1,3]
          # MR_result[temp,2] = presso_sd = presso_result$`Main MR results`[1,4]
          # MR_result[temp,3] = ifelse(presso_est-1.96*presso_sd<=beta_M&
          #                              presso_est+1.96*presso_sd>=beta_M,1,0)
          # MR_result[temp,4] = "MR-PRESSO"
          
        }else{
          
          
          
        }
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
      result$i_vec = rep(i,length(method))
      result$l_vec = rep(l,length(method))
      result$v_vec = rep(v,length(method))
      # result.list[[temp]] = result
      # temp = temp + 1
      # 
    #}
    
#   }
# }
  
 
#result = rbindlist(result.list)
save(result,file =paste0(cur.dir,"MR_result_chr22_beta_",i,"_rho_",l,"_ple_",v,".rdata"))

