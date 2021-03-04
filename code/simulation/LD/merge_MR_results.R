MR_result_list = list()
l = 3
  for(sub in 1:10){
    load(paste0(cur.dir,"MR_result_rho_",l,"_sub_",sub,".rdata"))
    MR_result_list[[sub]] = as.data.frame(MR_result)
  }

MR_result_com = rbindlist(MR_result_list)
MR_result_temp= MR_result_com %>% 
  group_by(method) %>% 
  summarize(mean_est = mean(est,na.rm =T),
            mean_cover = mean(cover,na.rm =T))
  mean(,na.rm=T)

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
  
  