i1 = 1
#i2 = as.numeric(args[[1]])
load("/data/zhangh24/MR_MA/result/simulation/prs/summary_gwas_our.rdata")
alpha_est_mat = result[[1]]
alpha_p_mat = result[[3]]
gamma_est_mat = result[[4]]




  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_Y.rdata"))
  
#genotype s2 is the genotype data for y
genotype_m_test = genotype_s2
dim(genotype_m_test)
n.rep = 1000
prs_m_mat = matrix(0,nrow(genotype_m_test),n.rep)
prs_y_mat = matrix(0,nrow(genotype_m_test),n.rep)
alpha_est = alpha_est_mat[,1]
# n.snp.cau = rep(0,length(pthres))
# n.snp.total = rep(0,length(pthres))

for(l in 1:n.rep){
  print(l)
  
    idx = c(1:5000)
  
    alpha_est_idx = alpha_est_mat[idx,l]
    gamma_est_idx = gamma_est_mat[idx,l]
    #find the corresponding SNPs within this genotype file  
    
    prs_m_mat[,l] = genotype_m_test%*%alpha_est_idx 
    prs_y_mat[,l] = genotype_m_test%*%gamma_est_idx
    
    
  }



prs_result = list(prs_m_mat,prs_y_mat)
save(prs_result,file =paste0("/data/zhangh24/MR_MA/result/simulation/prs/best_prs_result_twostage.rdata"))


# 
# 
# 
# (i1-1)*5000+1
# i1*5000
