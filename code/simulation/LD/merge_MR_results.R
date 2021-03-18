library(data.table)
library(dplyr)
MR_result_list = list()

cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
temp = 1
for(l in 1:3){
  for(sub in 1:10){
    load(paste0(cur.dir,"MR_result_rho_",l,"_sub_",sub,".rdata"))
    MR_result_list[[temp]] = data.frame(MR_result,l_vec = rep(l,nrow(MR_result)))
    temp  = temp+1
  }
  
}
 
MR_result_com = rbindlist(MR_result_list) %>% 
  mutate(est = as.numeric(est),
         sd = as.numeric(sd),
         cover = as.numeric(cover))
save(MR_result_com,file = paste0("/data/zhangh24/MR_MA/result/LD/MR_result_summary.rdata"))
  

  MR_result_com = as.data.frame(MR_result_com)
  idx <- which(MR_result_com$method=="IVW")
  mean(as.numeric(MR_result_com[idx,1]),na.rm=T)
  idx <- which(MR_result_com$method=="MR-Weight")
  mean(as.numeric(MR_result_com[idx,1]),na.rm=T)
  idx <- which(MR_result_com$method=="MR-Egger")
  mean(as.numeric(MR_result_com[idx,1]),na.rm=T)
  idx <- which(MR_result_com$method=="Raps")
  mean(as.numeric(MR_result_com[idx,1]),na.rm=T)
  idx <- which(MR_result_com$method=="MR-Median")
  mean(as.numeric(MR_result_com[idx,1]),na.rm=T)
  
  
  idx <- which(MR_result_com$method=="IVW")
  mean(as.numeric(MR_result_com[idx,3]),na.rm=T)
  idx <- which(MR_result_com$method=="MR-Weight")
  mean(as.numeric(MR_result_com[idx,3]),na.rm=T)
  idx <- which(MR_result_com$method=="MR-Egger")
  mean(as.numeric(MR_result_com[idx,3]),na.rm=T)
  idx <- which(MR_result_com$method=="Raps")
  mean(as.numeric(MR_result_com[idx,3]),na.rm=T)
  idx <- which(MR_result_com$method=="MR-Median")
  mean(as.numeric(MR_result_com[idx,3]),na.rm=T)
  
  