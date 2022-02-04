#p-value threshold 5E-04 and 1E-03 LD matrix uninvertable
#pthres = c(3E-06,1E-05,5E-05,1E-04,5E-04,1E-03)

cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
result.list = list()
r2_vec = c(0.001,0.2,0.4,0.6)
temp = 1
for(r_ind in 1:4){
  for(l in 1:3){
    load(paste0(cur.dir,"WMR_result_chr22_rho_",l,"_ple_",v,"r_ind",r_ind,".rdata"))
    result$method = paste0("WMR (LD r2 =",r2_vec[r_ind],")")
    result$l_vec = l
    result$v_vec = 1
    result.list[[temp]] = result  
    temp = temp + 1
  }
  
}
library(data.table)
result = rbindlist(result.list)
save(result,file = paste0(cur.dir,"WMR_result_chr22_LDgrid.rdata"))

