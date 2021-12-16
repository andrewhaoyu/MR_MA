library(data.table)
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
result.list = list()
temp =1
for(i in 1:2)
{
  for(l in 1:3){
    for(v in 1:3){
      load(paste0(cur.dir,"MR_result_chr22_beta_",i,"_rho_",l,"_ple_",v,".rdata"))
      result.list[[temp]] = result
      temp = temp + 1
    }
  }
}
result = rbindlist(result.list)
library(dplyr)
result %>% filter(
                    i_vec==1&
                    v_vec==1&
                      l_vec==3)
save(result,file = paste0(cur.dir,"MR_result_chr22.rdata"))
