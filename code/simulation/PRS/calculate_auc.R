#get the auc under different pthres
pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)
load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata"))
n.train = 100000
genotype_m_test = genotype_s[(n.train+1):nrow(genotype_s),]
prs_mat_all = matrix(0,nrow(genotype_m_test),length(pthres))
load("/data/zhangh24/MR_MA/result/simulation/prs/M_mat.rdata")
M_test = M_mat[(n.train+1):nrow(genotype_s),1]
for(i1 in 1:30){
  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/prs_mat_",i1,".rdata"))
  prs_mat_all = prs_mat_all+prs_mat
}

r2 = rep(0,length(pthres))
for(l in 1:length(pthres)){
  model =lm(M_test~prs_mat_all[,l])
  r2[l] = summary(model)$r.square
}