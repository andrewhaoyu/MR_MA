#Two different settings for reuslt name 
# load(paste0("./result/simulation/LD_simulation_test/result_np",i1,"_",i2,"_",i3,".rdata")) represent N 2000 p 1000 p_threshold 1E-5
# load(paste0("./result/simulation/LD_simulation_test/result",i1,"_",i2,"_",i3,".rdata")) represent N 60000 p 500 p_threshold 5E-08
setwd("/data/zhangh24/MR_MA/")
library(data.table)
beta_vec = c(0,0.2)
result.list = list()
temp = 1
pthres = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02)
n_pthres = length(pthres)
for(l in 1:3){
  for(v in 1:1){

    for(i1 in 1:length(pthres)){
      load(paste0("./result/simulation/LD_simulation_test/result_phtres_noLD",l,"_",v,"_",i1,".rdata"))
      #load(paste0("./result/simulation/LD_simulation_test/result_phtres_strongLD",l,"_",v,"_",i1,".rdata"))
      result.list[[temp]] = result  
      temp = temp + 1
    }
   
    
    
  }
  
}
result = rbindlist(result.list)
print(result)
save(result,file = paste0("./result/simulation/LD_simulation_test/wmr_simu_result_strongLD_pthres.rdata"))
