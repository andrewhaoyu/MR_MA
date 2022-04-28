args = commandArgs(trailingOnly = T)
#i for beta vector
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
#l is causal SNPs proportion: 0.05, 0.01, 0.001
#sub = as.numeric(args[[3]])

#r_ind for clumping grid
#r_ind 0.001, 0.2, 0.4, 0.6, 0.8, 1
#r_ind = as.numeric(args[[3]])
rep_ind = as.numeric(args[[3]])
#tau for pleotropic penalty
#tau as 0,1E-05,1E-04,1E-03,1E-02,1E-01
# tau_ind = as.numeric(args[[3]])
library(withr)
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/", install_github('qingyuanzhao/mr.raps'))
#install_github('qingyuanzhao/mr.raps')
#Goal: This analyses check the WMR performance under different penalty (tau_vec) with different LD

library(MASS)
library(MESS)
library(data.table)
library(dplyr)
library(Rfast)
library(mr.raps)
library(mr.divw)
# startend <- function(num,size,ind){
#   split.all <- split(1:num,cut(1:num,size))
#   temp <- split.all[[ind]]
#   start <- temp[1]
#   end <- temp[length(temp)]
#   return(c(start,end))
# }
# size = 25
# num = 100
# startend_result = startend(num,size,rep_ind)
# 
# start = startend_result[1]
# end = startend_result[2]

start = 1
end = 1

cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
setwd("/data/zhangh24/MR_MA/")
# l = 3
j = 22
#load LD score
temp = 1
result.list = list()
v =1
# tau_vec = c(0,1E-05,1E-04,1E-03,1E-02,1E-01)
# for(tau_ind in 1:length(tau_vec)){
#   tau = tau_vec[tau_ind]
# for(l in 1:3){
#   for(r_ind in 1:4){


#for(i in 1:2){
beta_vec = c(0,0.2)
beta = beta_vec[i]
#   for(sub in 1:10){
setwd("/data/zhangh24/MR_MA/")
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
j =22
sum.data.y = as.data.frame(fread(paste0(cur.dir,"y_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v)))
sum.data.m = as.data.frame(fread(paste0(cur.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v)))
sum.data.m2 = as.data.frame(fread(paste0(cur.dir,"m2_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v)))
#sum.data.m = left_join(sum.data.m,ldscore,by = c("ID"="SNP"))
#sum.data.m = inner_join(sum.data.m,ldscore,by = c("ID"="SNP"))
n.snp = nrow(sum.data.m)
n.rep = end-start+1
pthres = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02,1E-01,1)
r2_vec = c(0.001,0.2,0.4,0.6,0.8,1)
#pthres = c(1E-01,1)
for(r_ind in 1:1){
  for(i1 in 1:length(pthres)){
    beta_est = rep(0,n.rep)
    beta_cover = rep(0,n.rep)
    beta_se = rep(0,n.rep)
    
    i_rep = rep_ind
    #for(i_rep in  1:1){
    if(r_ind==6){
      #r_ind ==6 means no clumping at all
      LD.snp = sum.data.m[,"ID",drop=F]
      colnames(LD.snp) = "SNP"
      
    }else{
      LD.snp = as.data.frame(fread( paste0(cur.dir,"LD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep,"_rvec_",r_ind,".clumped")))  
    }
    
    
    
    sum.data.match.m = inner_join(LD.snp,sum.data.m,by = c("SNP"="ID"))
    matched.snp = sum.data.match.m[,"SNP",drop=F]
    p = sum.data.match.m[,(6+3*i_rep)]
    
    
    #idx = which(p<=pthres[i1])
    #idx = c(1,3,5)
    idx = which(p<=pthres[i1])
    select.snp = data.frame(SNP=sum.data.match.m$SNP[idx])
    #idx = 1
    #if(length(idx)>3){
    sum.data.match.y = left_join(matched.snp,sum.data.y,by=c("SNP"="ID"))
    sum.data.match.m2 = left_join(matched.snp,sum.data.m2,by=c("SNP"="ID"))
    Gamma = sum.data.match.y[,(6+3*i_rep-2)]
    se_Gamma = as.numeric(sum.data.match.y[,(6+3*i_rep-1)])
    alpha = as.numeric(sum.data.match.m2[,(6+3*i_rep-2)])
    se_alpha = as.numeric(sum.data.match.m2[,(6+3*i_rep-1)])
    SNP.select = matched.snp$SNP[idx]
    idx.match = match(SNP.select,sum.data.m$ID)
    
    alpha_select =alpha[idx]
    se_alpha_select = se_alpha[idx]
    Gamma_select = Gamma[idx]
    se_Gamma_select = se_Gamma[idx]
    
    data = data.frame(beta.exposure = alpha_select,
                      beta.outcome = Gamma_select,
                      se.exposure = se_alpha_select,
                      se.outcome = se_Gamma_select)
    
    raps_result = mr.raps(data,diagnostics = F)
    beta_est_Raps = raps_result$beta.hat
    beta_cover_Raps = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                                      raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
    beta_se_Raps = raps_result$beta.se

    result_raps = data.frame(beta_est = beta_est_Raps,
                             cover = beta_cover_Raps,
                             beta_se = beta_se_Raps,
                        i_vec = i,
                        p_vec = pthres[i1],
                        r_vec = r2_vec[r_ind],
                        method = "Raps")
    
    divw_result = mr.divw(alpha_select,
            Gamma_select,
            se_alpha_select,se_Gamma_select)
    
    beta_est_divw = divw_result$beta.hat
    beta_cover_divw = ifelse(divw_result$beta.hat-1.96*divw_result$beta.se<=beta&
                              divw_result$beta.hat+1.96*divw_result$beta.se>=beta,1,0)
    beta_se_divw =divw_result$beta.se
    
    result_divw = data.frame(beta_est = beta_est_divw,
                             cover = beta_cover_divw,
                             beta_se = beta_se_divw,
                             i_vec = i,
                             p_vec = pthres[i1],
                             r_vec = r2_vec[r_ind],
                             method = "divw")
    result = rbind(result_raps,result_divw)
    result.list[[temp]] = result
    temp = temp + 1
  }
  
  
  
  
}

#}


#}

#}


#   }
# }
result = rbindlist(result.list)
save(result,file = paste0(cur.dir,"other_method_result_chr22_beta_",i,"_rho_",l,"_rep_ind_",rep_ind,".rdata"))
# }
