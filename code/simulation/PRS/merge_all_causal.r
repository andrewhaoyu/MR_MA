#merge and plot mr test results
filedir = "/data/zhangh24/MR_MA/result/simulation/prs/"
#files = dir(filedir,"beta_test_result_500k_",full.names = T)
files = dir(filedir,"beta_test_result_rho_",full.names = T)

total = 0
for(i1 in 1:1000){
  #file  = paste0("/data/zhangh24/MR_MA/result/simulation/prs//beta_test_result_500k_",i1)
  file  = paste0("/data/zhangh24/MR_MA/result/simulation/prs//beta_test_result_rho_",i1)
  if(file %in% files){
    load(file)  
    temp = nrow(result[[1]])
    total = total+temp  
  }
  
}

n.rep = total
n.snp.vec = 3
beta_est_result = matrix(0,n.rep,n.snp.vec)
beta_est_result_cover = matrix(0,n.rep,n.snp.vec)
beta_est_result_summary = matrix(0,n.rep,n.snp.vec)
beta_est_result_summary_cover = matrix(0,n.rep,n.snp.vec)
beta_est_result_inner = matrix(0,n.rep,n.snp.vec)
beta_est_result_inner_cover = matrix(0,n.rep,n.snp.vec)
beta_est_result_inner_summary = matrix(0,n.rep,n.snp.vec)
beta_est_result_inner_cover_summary = matrix(0,n.rep,n.snp.vec)
beta_est_IVW = matrix(0,n.rep,n.snp.vec)
beta_est_IVW_cover = matrix(0,n.rep,n.snp.vec)
beta_est_IVW_inner = matrix(0,n.rep,n.snp.vec)
beta_est_IVW_inner_cover = matrix(0,n.rep,n.snp.vec)
total = 0
filedir = "/data/zhangh24/MR_MA/result/simulation/prs/"

for(i1 in 1:1000){
  file  = paste0("/data/zhangh24/MR_MA/result/simulation/prs//beta_test_result_rho_",i1)
  if(file %in% files){
    load(file)  
    temp = nrow(result[[1]])
    beta_est_result[total+(1:temp),] = result[[1]]
    beta_est_result_cover[total+(1:temp),] = result[[2]]
    beta_est_result_summary[total+(1:temp),] = result[[3]]
    beta_est_result_summary_cover[total+(1:temp),] = result[[4]]
    beta_est_result_inner[total+(1:temp),] = result[[5]]
    beta_est_result_inner_cover[total+(1:temp),] = result[[6]]
    beta_est_result_inner_summary[total+(1:temp),] = result[[7]]
    beta_est_result_inner_cover_summary[total+(1:temp),] = result[[8]]
    beta_est_IVW[total+(1:temp),] = result[[9]]
    beta_est_IVW_cover[total+(1:temp),] = result[[10]]
    beta_est_IVW_inner[total+(1:temp),] = result[[11]]
    beta_est_IVW_inner_cover[total+(1:temp),] = result[[12]]
    
    total = total+temp  
  }
  
}
Getestimate <- function(est_mat){
  return(paste0(sprintf("%.3f",round(colMeans(est_mat),3))," [",
                sprintf("%.3f",round(apply(est_mat,2,function(x){quantile(x,0.025)}),3)),
                ", ",
                sprintf("%.3f",round(apply(est_mat,2,function(x){quantile(x,0.975)}),3)),
                "]"))
  
  
}


method = rep(c("Best PRS","Best PRS (summary)","IVW (P<5E-08)"),2)
sample = c(rep("Two sample",3),rep("One sample",3))
result <- rbind(Getestimate(beta_est_result),
                Getestimate(beta_est_result_summary),
                Getestimate(beta_est_IVW),
                Getestimate(beta_est_result_inner),
                Getestimate(beta_est_result_inner_summary),
                Getestimate(beta_est_IVW_inner))

final.result <- cbind(sample,method,result)
write.csv(final.result, file = "/data/zhangh24/MR_MA/result/simulation/prs/prs_mr_result_rho0.3.csv",row.names = F)





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
s