cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
library(dplyr)
library(bc2)
library(data.table)
j = 22
num = 10


for(l in 1:3){
  i1 = 1
  sum.temp <- fread(paste0(cur.dir,"y_summary_chr_",j,"_rho_",l,"_sub_",i1))
  sum.infor <- sum.temp[,1:6]
  effect.list = list()
  for(i1 in 1:num){
    sum.temp <- fread(paste0(cur.dir,"y_summary_chr_",j,"_rho_",l,"_sub_",i1))
    effect.list[[i1]] = sum.temp[,7:ncol(sum.temp)]
  }
  effect <- bind_cols(effect.list)
  sum.data = cbind(sum.infor,effect)
  write.table(sum.data,file = paste0(cur.dir,"y_summary_chr_",j,"_rho_",l))
  effect.list = list()
  for(i1 in 1:num){
    sum.temp <- fread(paste0(cur.dir,"m_summary_chr_",j,"_rho_",l,"_sub_",i1))
    effect.list[[i1]] = sum.temp[,7:ncol(sum.temp)]
  }
  effect <- bind_cols(effect.list)
  sum.data = cbind(sum.infor,effect)
  write.table(sum.data,file = paste0(cur.dir,"m_summary_chr_",j,"_rho_",l))
}
system(paste0("rm ",cur.dir,"*_sub_*"))

# 
# for(l in 1:3){
#   sum.data = fread(paste0(cur.dir,"y_summary_chr_",j,"_rho_",l))
#   sum.data = sum.data[,-1]
#   write.table(sum.data,file = paste0(cur.dir,"y_summary_chr_",j,"_rho_",l),row.names = F,col.names = T,quote=F)
#   
#   sum.data = fread(paste0(cur.dir,"m_summary_chr_",j,"_rho_",l))
#   sum.data = sum.data[,-1]
#   
#   write.table(sum.data,file = paste0(cur.dir,"m_summary_chr_",j,"_rho_",l),row.names = F,col.names = T,quote=F)
# }
