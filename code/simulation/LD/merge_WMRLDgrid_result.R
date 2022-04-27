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
      result.list[[temp]] = result    
      temp = temp + 1
    }
  }
}
library(data.table)
result = rbindlist(result.list)

beta = beta_vec[i]
result.sub = result %>% 
  filter(i_vec==1&
           p_vec==5e-08&
           r_vec==0.001&
           l_vec == 1)
beta_est = result.sub$beta_est
se_est = result.sub$beta_se
result = data.frame(
  method = method,
  bias = apply(beta_est,2,function(x){mean(x,na.rm=T)})-beta,
  em_se = apply(beta_est,2,function(x){sd(x,na.rm=T)}),
  es_se = apply(se.result,2,function(x){mean(x,na.rm=T)}),
  cover = apply(cover.result,2,function(x){mean(x,na.rm=T)}),
  rmse = apply(mean.result,2,function(x){sqrt(mean((x-beta)^2,na.rm=T))})
)


beta_vec = c(0,0.2)
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
