args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)
load("/data/zhangh24/MR_MA/result/simulation/prs/summary_gwas_our.rdata")
alpha_est_mat = result[[1]]
alpha_p_mat = result[[3]]
gamma_est_mat = result[[4]]



if(i1==1){
  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata"))
  
}else{
  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/noncau_genotype_M_",i1,".rdata"))
  
}
n.train = 100000
genotype_m_test = genotype_s[(n.train+1):nrow(genotype_s),]
prs_mat = matrix(0,nrow(genotype_m_test),length(pthres))
alpha_est = alpha_est_mat[,1]
n.snp.cau = rep(0,length(pthres))
n.snp.total = rep(0,length(pthres))
for(l in 1:length(pthres)){
  
  idx <- which(alpha_p_mat[,1]<=pthres[l])
  n.snp.cau[l] = length(which(idx<=5000))
  n.snp.total[l] = length(idx)
  idx <- idx[idx>=((i1-1)*5000+1)&idx<=i1*5000]
  if(length(idx)>0){
    alpha_est_idx = alpha_est[idx]
    #find the corresponding SNPs within this genotype file  
    idx.scale = idx-(i1-1)*5000    
    prs_mat[,l] = genotype_m_test[,idx.scale,drop=F]%*%alpha_est_idx
  }

  
}

save(prs_mat,file =paste0("/data/zhangh24/MR_MA/result/simulation/prs/prs_mat_",i1,".rdata"))
# 
# 
# 
# (i1-1)*5000+1
# i1*5000
