args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])

pthres <- 1E-03
load("/data/zhangh24/MR_MA/result/simulation/prs/summary_gwas_our.rdata")
alpha_est_mat = result[[1]]
alpha_p_mat = result[[3]]
gamma_est_mat = result[[4]]



if(i1==1){
  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_Y.rdata"))
  
}else{
  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/noncau_genotype_Y_",i1,".rdata"))
  
}
#genotype s2 is the genotype data for y
genotype_m_test = genotype_s2
n.rep = 1000
prs_m_mat = matrix(0,nrow(genotype_m_test),n.rep)
prs_y_mat = matrix(0,nrow(genotype_m_test),n.rep)
alpha_est = alpha_est_mat[,1]
# n.snp.cau = rep(0,length(pthres))
# n.snp.total = rep(0,length(pthres))

  for(l in 1:n.rep){
    print(l)
    idx <- which(alpha_p_mat[,l]<=pthres)
    
    idx <- idx[idx>=((i1-1)*5000+1)&idx<=i1*5000]
    if(length(idx)>0){
      alpha_est_idx = alpha_est_mat[idx,l]
      gamma_est_idx = gamma_est_mat[idx,l]
      #find the corresponding SNPs within this genotype file  
      idx.scale = idx-(i1-1)*5000    
      prs_m_mat[,l] = genotype_m_test[,idx.scale,drop=F]%*%alpha_est_idx
      prs_y_mat[,l] = genotype_m_test[,idx.scale,drop=F]%*%gamma_est_idx
      
      
    }
  }
  

prs_result = list(prs_m_mat,prs_y_mat)
save(prs_result,file =paste0("/data/zhangh24/MR_MA/result/simulation/prs/prs_result_twostage_",i1,".rdata"))
# 
# 
# 
# (i1-1)*5000+1
# i1*5000
