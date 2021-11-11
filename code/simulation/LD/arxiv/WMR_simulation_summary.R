#Two different settings for reuslt name 
# load(paste0("./result/simulation/LD_simulation_test/result_np",i1,"_",i2,"_",i3,".rdata")) represent N 2000 p 1000 p_threshold 1E-5
# load(paste0("./result/simulation/LD_simulation_test/result",i1,"_",i2,"_",i3,".rdata")) represent N 60000 p 500 p_threshold 5E-08
setwd("/data/zhangh24/MR_MA/")
library(data.table)
beta_vec = c(1,0.5,0)
result.list = list()
temp = 1
for(i1 in 1:3){
  for(i2 in 1:3){
    mean.list = list()
    se.list = list()
    cover.list = list()
    for(i3 in 1:40){
      load(paste0("./result/simulation/LD_simulation_test/result_np",i1,"_",i2,"_",i3,".rdata"))
      mean.list[[i3]] = result[[1]]
      se.list[[i3]] = result[[2]]
      cover.list[[i3]] = result[[3]]
    }
    beta=beta_vec[i1]
    est = rbindlist(mean.list)
    se_est = rbindlist(se.list)
    cover = rbindlist(cover.list)
    result = data.frame(
      method = c("WMR","IVW","MR-Egger","MR-median","MRRAPs"),
      bias = apply(est,2,mean)-beta,
      em_se = apply(est,2,sd),
      es_se = apply(se_est,2,mean),
      cover = apply(cover,2,mean),
      rmse = apply(est,2,function(x){sqrt(mean((x-beta)^2))})
    )
    result$i1_vec = rep(i1,5)
    result$i2_vec = rep(i2,5)
    result.list[[temp]] = result
    temp = temp + 1
  }
  
}
result = rbindlist(result.list)
save(result,file = paste0("./result/simulation/LD_simulation_test/wmr_simu_result_com_np.rdata"))
