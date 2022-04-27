#p-value threshold 5E-04 and 1E-03 LD matrix uninvertable
pthres = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02,1E-01,1)
r2_vec = c(0.001,0.2,0.4,0.6,0.8,1)
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
result.list = list()
temp = 1
for(i in 1:2){
  for(l in 1:3){
    for(rep_ind in 1:100){
      load(paste0(cur.dir,"WMR_result_chr22_beta_",i,"_rho_",l,"_rep_ind_",rep_ind,".rdata"))
      result$l_vec = rep(l,nrow(result))
      result$r2_vec = rep(r2_vec,length(pthres))
      result.list[[temp]] = result    
      temp = temp + 1
    }
  }
}
library(data.table)
result = rbindlist(result.list)
save(result,file = paste0(cur.dir,"WMR_result_chr22_pthres.rdata"))



pthres = c(3E-06,1E-05,5E-05,1E-04,5E-04,1E-03)

cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
result.list = list()
for(i1 in 1:6){
  load(paste0(cur.dir,"WMR_result_no_LD_chr22_i1",i1,".rdata"))
  result$method = paste0("WMR (",pthres[i1],")")
  result.list[[i1]] = result
}
library(data.table)
result = rbindlist(result.list)
save(result,file = paste0(cur.dir,"WMR_result_no_LD_chr22_pthres.rdata"))
