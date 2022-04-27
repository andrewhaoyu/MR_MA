#p-value threshold 5E-04 and 1E-03 LD matrix uninvertable



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
