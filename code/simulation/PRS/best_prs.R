
load("/data/zhangh24/MR_MA/result/simulation/prs/summary_gwas_our.rdata")
alpha_est_mat = result[[1]]
alpha_p_mat = result[[3]]
gamma_est_mat = result[[4]]




load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata"))
  

n.train = 100000
genotype_m_test = genotype_s[(n.train+1):nrow(genotype_s),]
n.rep = 1000
prs_m_mat = matrix(0,nrow(genotype_m_test),n.rep)
prs_y_mat = matrix(0,nrow(genotype_m_test),n.rep)

# n.snp.cau = rep(0,length(pthres))
# n.snp.total = rep(0,length(pthres))

best_prs_est = rep(0,n.rep)
best_prs_low = rep(0,n.rep)
best_prs_high = rep(0,n.rep)
library(ISwR)
for(l in 1:n.rep){
  print(l)
  #idx <- c(1:5000)
  
  
    alpha_est_idx = alpha_est_mat[idx,l]
    gamma_est_idx = gamma_est_mat[idx,l]
    prs_m_mat[,l] = genotype_m_test%*%alpha_est_idx
    prs_y_mat[,l] = genotype_m_test%*%gamma_est_idx
    model <- lm(prs_y_mat[,2]~prs_m_mat[,2])
    #coefficients(summary(model))[2,1]
    temp_ci = confint(model)
    best_prs_est[l] = coefficients(summary(model))[2,1]
    best_prs_low[l] = temp_ci[2,1]
    best_prs_high[l] = temp_ci[2,2]
}





prs_est = rep(0,n.rep)
prs_low_est = rep(0,n.rep)
prs_high_est = rep(0,n.rep)
for(l in 1:1000){
  print(l)
  #idx <- c(1:5000)
  
  
  alpha_est_idx = alpha_est_mat[idx,l]
  gamma_est_idx = gamma_est_mat[idx,l]
  prs_m_mat[,l] = genotype_m_test[,idx]%*%alpha_est_idx
  prs_y_mat[,l] = genotype_m_test[,idx]%*%gamma_est_idx
  model <- lm(prs_y_mat[,l]~prs_m_mat[,l])
  #coefficients(summary(model))[2,1]
  temp_ci = confint(model)
  prs_est[l] = coefficients(summary(model))[2,1]
  prs_low_est[l] = temp_ci[2,1]
  prs_high_est[l] = temp_ci[2,2]
}


prs_result = list(prs_m_mat,prs_y_mat)
save(prs_result,file =paste0("/data/zhangh24/MR_MA/result/simulation/prs/best_prs_",i1,".rdata"))
# 
# 
# 
# (i1-1)*5000+1
# i1*5000
