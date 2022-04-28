#p-value threshold 5E-04 and 1E-03 LD matrix uninvertable
#pthres = c(3E-06,1E-05,5E-05,1E-04,5E-04,1E-03)
library(tidyverse)
pthres = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02,1E-01,1)
r2_vec = c(0.001,0.2,0.4,0.6,0.8,1)
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
result.list = list()
temp = 1
for(i in 1:2){
  for(l in 1:3){
    for(rep_ind in 1:100){
      load(paste0(cur.dir,"WMR_result_chr22_beta_",i,"_rho_",l,"_rep_ind_",rep_ind,".rdata"))
      result$l_vec = rep(l,nrow(result))
      result$rep_ind = rep(rep_ind,nrow(result))
      result.list[[temp]] = result    
      temp = temp + 1
    }
  }
}
library(data.table)
result = rbindlist(result.list)

result.table.list = list()
temp = 1

beta_vec = c(0,0.2)
for(i in 1:2){
  beta = beta_vec[i]
  for(l in 1:3){
    for(r in 1:length(r2_vec)){
      for(p in 1:length(pthres)){
        beta = beta_vec[i]
        result.sub = result %>% 
          filter(i_vec==i&
                   p_vec==pthres[p]&
                   r_vec==r2_vec[r]&
                   l_vec == l)
        beta_est = result.sub$beta_est
        se_est = result.sub$beta_se
        beta_cover = result.sub$beta_cover
        result.table = data.frame(
          method = method,
          bias = mean(beta_est)-beta,
          em_se = sd(beta_est),
          es_se = mean(beta_est),
          cover = mean(beta_cover),
          rmse = sqrt(mean((beta_est-beta)^2)),
          i_vec = i,
          p_vec = pthres[p],
          r_vec = r2_vec[r],
          l_vec = l
        )
        result.table.list[[temp]] = result.table
        temp = temp + 1
      }      
    }
  
  }
  }
sum.result = rbindlist(result.table.list)

result = sum.result
save(result,file = paste0(cur.dir,"WMR_result_chr22_grid.rdata"))

# mean.result = data.frame(
#   beta_est
# )
# colnames(mean.result) = method
# 
# se.result = data.frame(
#   beta_se
# )
# colnames(se.result) = method
# 
# cover.result = data.frame(
#   beta_cover
# )
# colnames(cover.result) = method
# 

# print(result)
# result$i_vec = rep(i,length(method))



save(result,file = paste0(cur.dir,"WMR_result_chr22_pthres.rdata"))