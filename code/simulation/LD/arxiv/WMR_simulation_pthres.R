args = commandArgs(trailingOnly = T)
l = as.numeric(args[[1]])
v = as.numeric(args[[2]])
i1 = as.numeric(args[[3]])

library(withr)
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/3.6/", install_github('qingyuanzhao/mr.raps'))
#install_github('qingyuanzhao/mr.raps')

library(MASS)
library(MESS)
library(data.table)
library(dplyr)
library(Rfast)


cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
setwd("/data/zhangh24/MR_MA/")
# l = 3
j = 22
#load LD score
source("./code/simulation/functions/MR_function_temp.R")
ldscore = fread("/data/zhangh24/MR_MA/data/eur_w_ld_chr/22.l2.ldscore.gz")
ldscore = ldscore %>% 
  dplyr::select(SNP,L2)
#load SNP id match file
load("/data/zhangh24/MR_MA/result/LD/chr22_snp_infor.rdata")
ldscore = left_join(snp.infor.subset,ldscore,
                    by=c("rs_id"="SNP"))
#load LD matrix
load(paste0(cur.dir,"chr_",j,"_LDmat.rdata"))
temp = 1
result.list = list()



for(i in 1:2){
  beta_vec = c(0,0.2)
  beta = beta_vec[i]
  #   for(sub in 1:10){
  setwd("/data/zhangh24/MR_MA/")
  cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
  j =22
  sum.data.y = as.data.frame(fread(paste0(cur.dir,"y_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v)))
  sum.data.m = as.data.frame(fread(paste0(cur.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v))) 
  sum.data.m = left_join(sum.data.m,ldscore,by = c("ID"="SNP"))
  n.snp = nrow(sum.data.m)
  n.rep = 100
  
  beta_est = rep(0,n.rep)
  beta_cover = rep(0,n.rep)
  beta_se = rep(0,n.rep)
  pthres = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02)
  for(i_rep in  1:n.rep){
    
    LD.snp = as.data.frame(fread( paste0(cur.dir,"LD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep,".clumped")))
    sum.data.match.m = left_join(LD.snp,sum.data.m,by = c("SNP"="ID"))
    p = sum.data.match.m[,(6+3*i_rep)]
    
    
    #idx = which(p<=pthres[i1])
    #idx = c(1,3,5)
    idx = which(p<=pthres[i1])
    #idx = 1
    #if(length(idx)>3){
    sum.data.match.y = left_join(LD.snp,sum.data.y,by=c("SNP"="ID"))
    Gamma = sum.data.match.y[,(6+3*i_rep-2)]
    se_Gamma = as.numeric(sum.data.match.y[,(6+3*i_rep-1)])
    alpha = as.numeric(sum.data.match.m[,(6+3*i_rep-2)])
    se_alpha = as.numeric(sum.data.match.m[,(6+3*i_rep-1)])
    MAF = sum.data.m[,"MAF"]
    SNP.select = sum.data.match.m$SNP[idx]
    idx.match = match(SNP.select,sum.data.m$ID)
    
    alpha_select =alpha[idx]
    se_alpha_select = se_alpha[idx]
    Gamma_select = Gamma[idx]
    se_Gamma_select = se_Gamma[idx]
    ld_score_select = ldscore$L2[idx.match]
    R = as.matrix(corr0[idx.match,idx.match])
    MAF_select = MAF[idx.match]
    
    
    MR_result <- WMRFun(Gamma_select,se_Gamma_select,
                        alpha_select,se_alpha_select,
                        ld_score_select,R,MAF_select)
    # MRWeight(Gamma = sumGamma,
    #                     var_Gamma = var_Gamma,
    #                     alpha = sumalpha,
    #                     var_alpha = var_alpha,
    #                     R = R)
    beta_est[i_rep] = MR_result[1]
    beta_cover[i_rep] = ifelse(MR_result[3]<=beta&MR_result[4]>=beta,1,0)
    beta_se[i_rep] = MR_result[2]
    
  }
  
  method = paste0("WMR (p<",pthres[i1],")")
  
  mean.result = data.frame(
    beta_est
  )
  colnames(mean.result) = method
  
  se.result = data.frame(
    beta_se
  )
  colnames(se.result) = method
  
  cover.result = data.frame(
    beta_cover
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
  result.list[[temp]] = result
  temp = temp + 1
  
}

result = rbindlist(result.list)
save(result,file = paste0("./result/simulation/LD_simulation_test/result_phtres_noLD",l,"_",v,"_",i1,".rdata"))
# }
