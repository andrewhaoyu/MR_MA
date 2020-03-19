#Merge the ratio estiamtes results
setwd("/data/zhangh24/MR_MA/")

n_vec = c(15000,75000,150000)
alpha_vec = c(0.0,0.01,0.03,0.05)
beta_vec = c(0,0.2,0.5)
times = 1000*100
replicates = 100
result_final = list()
temp.idx = 1
for(i4 in 1:3){
for(i1 in 1:3){
  for(i2 in 1:4){
    
     
      Gamma_est = rep(0,times)
      Gamma_var = rep(0,times)
      gamma_est = rep(0,times)
      gamma_var = rep(0,times)
      ratio_est = rep(0,times)
      ratio_var_standard = rep(0,times)
      ratio_var_alpha = rep(0,times)
      cover_ratio = rep(0,times)
      cover_alpha = rep(0,times)
      ci_low_ratio = rep(0,times)
      ci_high_ratio = rep(0,times)
      ci_low_alpha = rep(0,times)
      ci_high_alpha = rep(0,times)
      total = 0
      for(i3 in 1:100){
        load(paste0("./result/simulation/ratio_estimate/ratio_estimate_delta_",i1,"_",i2,"_",i3,"_",i4,".Rdata"))
        temp = length(result[[1]])
        Gamma_est[total+(1:temp)] = result[[1]]
        Gamma_var[total+(1:temp)] = result[[2]]
        gamma_est [total+(1:temp)] = result[[3]]
        gamma_var[total+(1:temp)] = result[[4]]
        ratio_est[total+(1:temp)] = result[[5]]
        ratio_var_standard[total+(1:temp)] = result[[6]]
        ratio_var_alpha[total+(1:temp)] =result[[7]]
        cover_ratio[total+(1:temp)] = result[[8]]
        cover_alpha[total+(1:temp)] = result[[9]]
        ci_low_ratio[total+(1:temp)] = result[[10]]
        ci_high_ratio[total+(1:temp)] = result[[11]]
        ci_low_alpha[total+(1:temp)] = result[[12]]
        ci_high_alpha[total+(1:temp)] = result[[13]]
        
        total = total+temp
        
        
      }
      #two loops
      #first loop sample size
      #inner loop for alpha_G
      result_final[[temp.idx]] =  list(Gamma_est,
                                       Gamma_var,
                                       gamma_est,
                                       gamma_var,
                                       ratio_est,
                                       ratio_var_standard,
                                       ratio_var_alpha,
                                       cover_ratio,
                                       cover_alpha,
                                       ci_low_ratio,
                                       ci_high_ratio,
                                       ci_low_alpha,
                                       ci_high_alpha)
      temp.idx = temp.idx+1
    }
    
    
  }
}
save(result_final,file = paste0("./result/simulation/ratio_estimate/ratio_estimate_merged_delta.Rdata"))



