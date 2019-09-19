#goal: merge the simulations datasets
times = 10
reps = 1000
TwoStage_est = rep(0,times*reps)
IVW_est = rep(0,times*reps)
IVWs_est = rep(0,times*reps)
cover_TwoStage_est = rep(0,times*reps)
cover_IVW_est = rep(0,times*reps)
cover_IVWs_est = rep(0,times*reps)
sigma_TwoStage = rep(0,times*reps)
sigma_IVW = rep(0,times*reps)
sigma_IVWs = rep(0,times*reps)
twostage.nsnps = rep(0,times*reps)
twostage.prop = rep(0,times*reps)
IVW.nsnps = rep(0,times*reps)
IVW.prop = rep(0,times*reps)

setwd("/spin1/users/zhangh24/MR_MA/")

total = 0
for(i1 in 1:reps){
  load(paste0("./result/simulation/simulation_",i1,".Rdata"))
  temp = length(result[[1]][[1]])
  TwoStage_est[total+(1:temp)] <- result[[1]][[1]]
  IVW_est[total+(1:temp)] <- result[[1]][[2]]
  IVWs_est[total+(1:temp)] <- result[[1]][[3]]
  cover_TwoStage_est[total+(1:temp)] <- result[[1]][[4]]
  cover_IVW_est[total+(1:temp)] <- result[[1]][[5]]
  cover_IVWs_est[total+(1:temp)] = result[[1]][[6]]
  sigma_TwoStage[total+(1:temp)] = result[[1]][[7]]
  sigma_IVW[total+(1:temp)] = result[[1]][[8]]
  sigma_IVWs[total+(1:temp)] = result[[1]][[9]]
  total = total+temp
}
IVWs_est_new = IVWs_est[order(IVWs_est)][(times*reps*0.05):(times*reps*0.95)]
beta_M = 0.1
mean(TwoStage_est)-beta_M
mean(IVW_est)-beta_M
mean(IVWs_est)-beta_M
mean(cover_TwoStage_est)
mean(cover_IVW_est)
mean(cover_IVWs_est,na.rm = T)
var(TwoStage_est)
var(IVW_est)
var(IVWs_est)
mean(sigma_TwoStage)
mean(sigma_IVW,na.rm = T)
mean(sigma_IVWs[order(IVWs_est)][(times*reps*0.05):(times*reps*0.95)])








total = 0
for(i1 in 1:reps){
  load(paste0("./result/simulation/simulation_",i1,".Rdata"))
  temp = length(result[[3]][[1]])
  TwoStage_est[total+(1:temp)] <- result[[3]][[1]]
  IVW_est[total+(1:temp)] <- result[[3]][[2]]
  IVWs_est[total+(1:temp)] <- result[[3]][[3]]
  cover_TwoStage_est[total+(1:temp)] <- result[[3]][[4]]
  cover_IVW_est[total+(1:temp)] <- result[[3]][[5]]
  cover_IVWs_est[total+(1:temp)] = result[[3]][[6]]
  sigma_TwoStage[total+(1:temp)] = result[[3]][[7]]
  sigma_IVW[total+(1:temp)] = result[[3]][[8]]
  sigma_IVWs[total+(1:temp)] = result[[3]][[9]]
  total = total+temp
}
beta_M = 0.1
mean(TwoStage_est)-beta_M
mean(IVW_est)-beta_M
mean(IVWs_est)-beta_M
IVWs_est_new = IVWs_est[order(IVWs_est)][(times*reps*0.05):(times*reps*0.95)]
mean(IVWs_est_new)
var(IVWs_est_new)
mean(cover_TwoStage_est)
mean(cover_IVW_est)
mean(cover_IVWs_est,na.rm = T)
var(TwoStage_est)
var(IVW_est)
var(IVWs_est)
mean(sigma_TwoStage)
mean(sigma_IVW,na.rm = T)
mean(sigma_IVWs)






total = 0
n_thres = 7
TwoStage_est_all <- matrix(0,times*reps,n_thres)
IVW_est_all <- matrix(0,times*reps,n_thres)
IVWs_est_all <- matrix(0,times*reps,n_thres)

for(i1 in 1:reps){
  load(paste0("./result/simulation/simulation_",i1,".Rdata"))
  temp = length(result[[2]][[1]])
  TwoStage_est[total+(1:temp)] <- result[[2]][[1]]
  IVW_est[total+(1:temp)] <- result[[2]][[2]]
  IVWs_est[total+(1:temp)] <- result[[2]][[3]]
  cover_TwoStage_est[total+(1:temp)] <- result[[2]][[4]]
  cover_IVW_est[total+(1:temp)] <- result[[2]][[5]]
  cover_IVWs_est[total+(1:temp)] = result[[2]][[6]]
  sigma_TwoStage[total+(1:temp)] = result[[2]][[7]]
  sigma_IVW[total+(1:temp)] = result[[2]][[8]]
  sigma_IVWs[total+(1:temp)] = result[[2]][[9]]
  twostage.nsnps[total+(1:temp)] = result[[2]][[10]]
  twostage.prop[total+(1:temp)] = result[[2]][[11]]
  IVW.nsnps[total+(1:temp)] = result[[2]][[12]]
  IVW.prop[total+(1:temp)] = result[[2]][[13]]
  TwoStage_est_all[total+(1:temp),] = result[[2]][[14]]
  IVW_est_all[total+(1:temp),] = result[[2]][[15]]
  IVWs_est_all[total+(1:temp),] = result[[2]][[16]]
  total = total+temp
}
beta_M = 0.1
mean(TwoStage_est)-beta_M
mean(IVW_est,na.rm = T)-beta_M
mean(IVWs_est,na.rm = T)-beta_M
mean(cover_TwoStage_est,na.rm = T)
mean(cover_IVW_est,na.rm = T)
mean(cover_IVWs_est,na.rm = T)
var(TwoStage_est,na.rm = T)
var(IVW_est,na.rm = T)
var(IVWs_est,na.rm = T)
mean(sigma_TwoStage,na.rm = T)
mean(sigma_IVW,na.rm = T)
mean(sigma_IVWs,na.rm = T)
mean(twostage.nsnps,na.rm = T)
mean(twostage.prop,na.rm = T)
mean(IVW.nsnps,na.rm = T)
mean(IVW.prop)








total = 0
for(i1 in 1:reps){
  load(paste0("./result/simulation/simulation_",i1,".Rdata"))
  temp = length(result[[4]][[1]])
  TwoStage_est[total+(1:temp)] <- result[[4]][[1]]
  IVW_est[total+(1:temp)] <- result[[4]][[2]]
  IVWs_est[total+(1:temp)] <- result[[4]][[3]]
  cover_TwoStage_est[total+(1:temp)] <- result[[4]][[4]]
  cover_IVW_est[total+(1:temp)] <- result[[4]][[5]]
  cover_IVWs_est[total+(1:temp)] = result[[4]][[6]]
  sigma_TwoStage[total+(1:temp)] = result[[4]][[7]]
  sigma_IVW[total+(1:temp)] = result[[4]][[8]]
  sigma_IVWs[total+(1:temp)] = result[[4]][[9]]
  twostage.nsnps[total+(1:temp)] = result[[4]][[10]]
  twostage.prop[total+(1:temp)] = result[[4]][[11]]
  IVW.nsnps[total+(1:temp)] = result[[4]][[12]]
  IVW.prop[total+(1:temp)] = result[[4]][[13]]
  
  total = total+temp
}
beta_M = 0.1
mean(TwoStage_est)-beta_M
mean(IVW_est,na.rm = T)-beta_M
mean(IVWs_est,na.rm = T)-beta_M
mean(cover_TwoStage_est,na.rm = T)
mean(cover_IVW_est,na.rm = T)
mean(cover_IVWs_est,na.rm = T)
var(TwoStage_est,na.rm = T)
var(IVW_est,na.rm = T)
var(IVWs_est,na.rm = T)
mean(sigma_TwoStage,na.rm = T)
mean(sigma_IVW,na.rm = T)
mean(sigma_IVWs,na.rm = T)
mean(twostage.nsnps,na.rm = T)
mean(twostage.prop,na.rm = T)
mean(IVW.nsnps,na.rm = T)
mean(IVW.prop,na.rm = T)

