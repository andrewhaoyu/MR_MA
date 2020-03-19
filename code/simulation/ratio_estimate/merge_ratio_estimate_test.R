#Merge the ratio estiamtes results
setwd("/data/zhangh24/MR_MA/")

n_vec = c(15000,75000,150000)
alpha_vec = c(0.0,0.01,0.03,0.05)
times = 1000*100
replicates = 100
result_final = list()
temp.idx = 1
for(i1 in 1:3){
  for(i2 in 1:4){
    
    Gamma_est = rep(0,times)
    Gamma_var = rep(0,times)
    gamma_est = rep(0,times)
    gamma_var = rep(0,times)
    ratio_est = rep(0,times)
    ratio_var = rep(0,times)
    ratio_p1 = rep(0,times)
    ratio_p2 = rep(0,times)
    p_weight = rep(0,times)
    p_delta = rep(0,times)
    # idx = which(ratio_est>=q_result[1]&
    #                ratio_est<=-6.24)
    # head(ci_low_exact[idx])
    total = 0
    for(i3 in 1:100){
      load(paste0("./result/simulation/ratio_estimate/p_estimate_",i1,"_",i2,"_",i3,".Rdata"))
      temp = length(result[[1]])
      Gamma_est[total+(1:temp)] = result[[1]]
      Gamma_var[total+(1:temp)] = result[[2]]
      gamma_est [total+(1:temp)] = result[[3]]
      gamma_var[total+(1:temp)] = result[[4]]
      ratio_est[total+(1:temp)] = result[[5]]
      ratio_var[total+(1:temp)] = result[[6]]
      ratio_p1[total+(1:temp)] =result[[7]]
      ratio_p2[total+(1:temp)] = result[[8]]
      p_weight[total+(1:temp)] = result[[9]]
      p_delta[total+(1:temp)] = result[[10]]
      total = total+temp
      
     
    }
    #two loops
    #first loop sample size
    #inner loop for alpha_G
    result_final[[temp.idx]] = list(Gamma_est,
                                     Gamma_var,
                                     gamma_est,
                                     gamma_var,
                                     ratio_est,
                                     ratio_var,
                                    ratio_p1,
                                    ratio_p2,
                                    p_weight,
                                    p_delta)
    temp.idx = temp.idx+1
  }
}
save(result_final,file = paste0("./result/simulation/ratio_estimate/ratio_estimate_merged_test.Rdata"))



