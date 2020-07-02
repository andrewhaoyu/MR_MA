#merge and plot mr test results
result= list(beta_est_result,
             beta_est_IVW,
             beta_est_MRLR,
             beta_est_result_inner,
             beta_est_IVW_inner,
             beta_est_MRLR_inner)
files = dir(filedir,"beta_test_result_",full.names = T)
n.rep = 2*length(files)
n.snp.vec = 3
beta_est_result = matrix(0,n.rep,n.snp.vec)
beta_est_IVW = matrix(0,n.rep,n.snp.vec)
beta_est_MRLR = matrix(0,n.rep,n.snp.vec)
beta_est_result_inner = matrix(0,n.rep,n.snp.vec)
beta_est_IVW_inner = matrix(0,n.rep,n.snp.vec)
beta_est_MRLR_inner = matrix(0,n.rep,n.snp.vec)
total = 0
filedir = "/data/zhangh24/MR_MA/result/simulation/prs/"

for(i1 in 1:1000){
  file  = paste0("/data/zhangh24/MR_MA/result/simulation/prs//beta_test_result_",i1)
  if(file %in% files){
    load(file)  
    temp = nrow(result[[1]])
    beta_est_result[total+(1:temp),] = result[[1]]
    beta_est_IVW[total+(1:temp),] = result[[2]]
    beta_est_MRLR[total+(1:temp),] = result[[3]]
    beta_est_result_inner[total+(1:temp),] = result[[4]]
    beta_est_IVW_inner[total+(1:temp),] = result[[5]]
    beta_est_MRLR_inner[total+(1:temp),] = result[[6]]
    total = total+temp  
  }
  
}

method = c("Best PRS","IVW (P<5E-08)","MR-weighted (P<5E-08)")
est_beta_est = colMeans(beta_est_result)
colMeans(beta_est_IVW)
colMeans(beta_est_MRLR)
colMeans(beta_est_result_inner)
colMeans(beta_est_IVW_inner)
colMeans(beta_est_MRLR_inner)
round(apply(beta_est_result,2,function(x){quantile(x,0.025)}),3)
round(apply(beta_est_result,2,function(x){quantile(x,0.975)}),3)

round(apply(beta_est_IVW,2,function(x){quantile(x,0.025)}),3)
round(apply(beta_est_IVW,2,function(x){quantile(x,0.975)}),3)
round(apply(beta_est_MRLR,2,function(x){quantile(x,0.025)}),3)
round(apply(beta_est_MRLR,2,function(x){quantile(x,0.975)}),3)
round(apply(beta_est_result_inner,2,function(x){quantile(x,0.025)}),3)
round(apply(beta_est_result_inner,2,function(x){quantile(x,0.975)}),3)
round(apply(beta_est_IVW_inner,2,function(x){quantile(x,0.025)}),3)
round(apply(beta_est_IVW_inner,2,function(x){quantile(x,0.975)}),3)
round(apply(beta_est_MRLR_inner,2,function(x){quantile(x,0.025)}),3)
round(apply(beta_est_MRLR_inner,2,function(x){quantile(x,0.975)}),3)
