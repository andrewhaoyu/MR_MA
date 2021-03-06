#goal: merge the simulations datasets
Genearte_result <- function(i2){
  times = 100
  reps = 100
  TwoStage_est = rep(0,times*reps)
  cover_TwoStage_est = rep(0,times*reps)
  sigma_TwoStage = rep(0,times*reps)
  sigma_y_TwoStage = rep(0,times*reps)
  IVW_est = rep(0,times*reps)
  cover_IVW_est = rep(0,times*reps)
  sigma_IVW = rep(0,times*reps)
  sigma_y_IVW = rep(0,times*reps)
  IVW_meta_est = rep(0,times*reps)
  cover_IVW_meta_est = rep(0,times*reps)
  sigma_IVW_meta = rep(0,times*reps)
  sigma_y_IVW_meta = rep(0,times*reps)
  IVWs_est = rep(0,times*reps)
  cover_IVWs_est = rep(0,times*reps)
  sigma_IVWs = rep(0,times*reps)
  setwd("/spin1/users/zhangh24/MR_MA/")
  total = 0
  for(i1 in 1:reps){
    load(paste0("./result/simulation/simulation_M_fixed",i1,"_",i2,".Rdata"))
    temp = length(result[[1]][[1]])
    TwoStage_est[total+(1:temp)] <- result[[1]][[1]]
    cover_TwoStage_est[total+(1:temp)] = result[[1]][[2]]
    sigma_TwoStage[total+(1:temp)] = result[[1]][[3]]
    sigma_y_TwoStage[total+(1:temp)] = result[[1]][[4]]
    IVW_est[total+(1:temp)] = result[[1]][[5]]
    cover_IVW_est[total+(1:temp)] = result[[1]][[6]]
    sigma_IVW[total+(1:temp)] = result[[1]][[7]]
    sigma_y_IVW[total+(1:temp)] = result[[1]][[8]]
    IVW_meta_est[total+(1:temp)] = result[[1]][[9]]
    cover_IVW_meta_est[total+(1:temp)] = result[[1]][[10]]
    sigma_IVW_meta[total+(1:temp)] = result[[1]][[11]]
    sigma_y_IVW_meta[total+(1:temp)] = result[[1]][[12]]
    IVWs_est[total+(1:temp)] = result[[1]][[13]]
    cover_IVWs_est[total+(1:temp)] = result[[1]][[14]]
    sigma_IVWs[total+(1:temp)] = result[[1]][[15]]
    total = total+temp
  }
  # IVWs_est_new = IVWs_est[order(IVWs_est)][(times*reps*0.05):(times*reps*0.95)]
  beta_M = 0.1
  
  bias <- rep(NA,4)
  bias[1] <- mean(TwoStage_est)-beta_M
  bias[2] <- mean(IVW_est)-beta_M
  bias[3] <- mean(IVW_meta_est)-beta_M
  bias[4] <- mean(IVWs_est)-beta_M
  EM_var <- rep(0,4)
  EM_var[1] <- var(TwoStage_est)
  EM_var[2] <- var(IVW_est)
  EM_var[3] <- var(IVW_meta_est)
  EM_var[4] <- var(IVWs_est)
  est_var <- rep(0,4)
  est_var[1] <- mean(sigma_TwoStage)
  est_var[2] <- mean(sigma_IVW)
  est_var[3] <- mean(sigma_IVW_meta)
  est_var[4] <- mean(sigma_IVWs)
  cover <- rep(0,4)
  cover[1] <- mean(cover_TwoStage_est)
  cover[2] <- mean(cover_IVW_est)
  cover[3] <- mean(cover_IVW_meta_est,na.rm=T)
  cover[4] <- mean(cover_IVWs_est,na.rm = T)
  
  est_var_y <- rep(NA,4)
  est_var_y[1] <- mean(sigma_y_TwoStage)
  est_var_y[2] <- mean(sigma_y_IVW)
  est_var_y[3] <- mean(sigma_y_IVW_meta)
  
  data1 <- data.frame(bias,EM_var,est_var,cover,est_var_y)
  
  
  times = 100
  reps = 100
  TwoStage_est = rep(0,times*reps)
  cover_TwoStage_est = rep(0,times*reps)
  sigma_TwoStage = rep(0,times*reps)
  sigma_y_TwoStage = rep(0,times*reps)
  IVW_est = rep(0,times*reps)
  cover_IVW_est = rep(0,times*reps)
  sigma_IVW = rep(0,times*reps)
  sigma_y_IVW = rep(0,times*reps)
  IVW_meta_est = rep(0,times*reps)
  cover_IVW_meta_est = rep(0,times*reps)
  sigma_IVW_meta = rep(0,times*reps)
  sigma_y_IVW_meta = rep(0,times*reps)
  IVWs_est = rep(0,times*reps)
  cover_IVWs_est = rep(0,times*reps)
  sigma_IVWs = rep(0,times*reps)
  setwd("/spin1/users/zhangh24/MR_MA/")
  total = 0
  for(i1 in 1:reps){
    load(paste0("./result/simulation/simulation_M_fixed",i1,"_",i2,".Rdata"))
    temp = length(result[[2]][[1]])
    TwoStage_est[total+(1:temp)] <- result[[2]][[1]]
    cover_TwoStage_est[total+(1:temp)] = result[[2]][[2]]
    sigma_TwoStage[total+(1:temp)] = result[[2]][[3]]
    sigma_y_TwoStage[total+(1:temp)] = result[[2]][[4]]
    IVW_est[total+(1:temp)] = result[[2]][[5]]
    cover_IVW_est[total+(1:temp)] = result[[2]][[6]]
    sigma_IVW[total+(1:temp)] = result[[2]][[7]]
    sigma_y_IVW[total+(1:temp)] = result[[2]][[8]]
    IVW_meta_est[total+(1:temp)] = result[[2]][[9]]
    cover_IVW_meta_est[total+(1:temp)] = result[[2]][[10]]
    sigma_IVW_meta[total+(1:temp)] = result[[2]][[11]]
    sigma_y_IVW_meta[total+(1:temp)] = result[[2]][[12]]
    IVWs_est[total+(1:temp)] = result[[2]][[13]]
    cover_IVWs_est[total+(1:temp)] = result[[2]][[14]]
    sigma_IVWs[total+(1:temp)] = result[[2]][[15]]
    total = total+temp
  }
  # IVWs_est_new = IVWs_est[order(IVWs_est)][(times*reps*0.05):(times*reps*0.95)]
  
  bias <- rep(NA,4)
  bias[1] <- mean(TwoStage_est)-beta_M
  bias[2] <- mean(IVW_est)-beta_M
  bias[3] <- mean(IVW_meta_est)-beta_M
  bias[4] <- mean(IVWs_est)-beta_M
  EM_var <- rep(0,4)
  EM_var[1] <- var(TwoStage_est)
  EM_var[2] <- var(IVW_est)
  EM_var[3] <- var(IVW_meta_est)
  EM_var[4] <- var(IVWs_est)
  est_var <- rep(0,4)
  est_var[1] <- mean(sigma_TwoStage)
  est_var[2] <- mean(sigma_IVW)
  est_var[3] <- mean(sigma_IVW_meta)
  est_var[4] <- mean(sigma_IVWs)
  cover <- rep(0,4)
  cover[1] <- mean(cover_TwoStage_est)
  cover[2] <- mean(cover_IVW_est)
  cover[3] <- mean(cover_IVW_meta_est,na.rm=T)
  cover[4] <- mean(cover_IVWs_est,na.rm = T)
  
  est_var_y <- rep(NA,4)
  est_var_y[1] <- mean(sigma_y_TwoStage)
  est_var_y[2] <- mean(sigma_y_IVW)
  est_var_y[3] <- mean(sigma_y_IVW_meta)
    
  data2 <- data.frame(bias,EM_var,est_var,cover,est_var_y)
 
  
  data <- rbind(data1,data2)
  return(data)
}
data <- Genearte_result(1)
data <- round(data,4)
write.csv(data,file = "./result/simulation/beta_fixed_g_0.01_result.csv")
data <- Genearte_result(2)
data <- round(data,5)
write.csv(data,file = "./result/simulation/beta_fixed_g_0.05_result.csv")

times = 100
reps = 100
TwoStage_est = rep(0,times*reps)
cover_TwoStage_est = rep(0,times*reps)
sigma_TwoStage = rep(0,times*reps)
sigma_y_TwoStage = rep(0,times*reps)
IVW_est = rep(0,times*reps)
cover_IVW_est = rep(0,times*reps)
sigma_IVW = rep(0,times*reps)
sigma_y_IVW = rep(0,times*reps)
IVW_meta_est = rep(0,times*reps)
cover_IVW_meta_est = rep(0,times*reps)
sigma_IVW_meta = rep(0,times*reps)
sigma_y_IVW_meta = rep(0,times*reps)
IVWs_est = rep(0,times*reps)
cover_IVWs_est = rep(0,times*reps)
sigma_IVWs = rep(0,times*reps)
setwd("/spin1/users/zhangh24/MR_MA/")
i2 = 2
total = 0
for(i1 in 1:reps){
  load(paste0("./result/simulation/simulation_",i1,"_",i2,".Rdata"))
  temp = length(result[[1]][[1]])
  TwoStage_est[total+(1:temp)] <- result[[1]][[1]]
  cover_TwoStage_est[total+(1:temp)] = result[[1]][[2]]
  sigma_TwoStage[total+(1:temp)] = result[[1]][[3]]
  sigma_y_TwoStage[total+(1:temp)] = result[[1]][[4]]
  IVW_est[total+(1:temp)] = result[[1]][[5]]
  cover_IVW_est[total+(1:temp)] = result[[1]][[6]]
  sigma_IVW[total+(1:temp)] = result[[1]][[7]]
  sigma_y_IVW[total+(1:temp)] = result[[1]][[8]]
  IVW_meta_est[total+(1:temp)] = result[[1]][[9]]
  cover_IVW_meta_est[total+(1:temp)] = result[[1]][[10]]
  sigma_IVW_meta[total+(1:temp)] = result[[1]][[11]]
  sigma_y_IVW_meta[total+(1:temp)] = result[[1]][[12]]
  IVWs_est[total+(1:temp)] = result[[1]][[13]]
  cover_IVWs_est[total+(1:temp)] = result[[1]][[14]]
  sigma_IVWs[total+(1:temp)] = result[[1]][[15]]
  total = total+temp
}
# IVWs_est_new = IVWs_est[order(IVWs_est)][(times*reps*0.05):(times*reps*0.95)]
beta_M = 0.1

bias <- rep(NA,4)
bias[1] <- mean(TwoStage_est)-beta_M
bias[2] <- mean(IVW_est)-beta_M
bias[3] <- mean(IVW_meta_est)-beta_M
bias[4] <- mean(IVWs_est)-beta_M
EM_var <- rep(0,4)
EM_var[1] <- var(TwoStage_est)
EM_var[2] <- var(IVW_est)
EM_var[3] <- var(IVW_meta_est)
EM_var[4] <- var(IVWs_est)
est_var <- rep(0,4)
est_var[1] <- mean(sigma_TwoStage)
est_var[2] <- mean(sigma_IVW)
est_var[3] <- mean(sigma_IVW_meta)
est_var[4] <- mean(sigma_IVWs)
cover <- rep(0,4)
cover[1] <- mean(cover_TwoStage_est)
cover[2] <- mean(cover_IVW_est)
cover[3] <- mean(cover_IVW_meta_est,na.rm=T)
cover[4] <- mean(cover_IVWs_est,na.rm = T)

est_var_y <- rep(NA,4)
est_var_y[1] <- mean(sigma_y_TwoStage)
est_var_y[2] <- mean(sigma_y_IVW)
est_var_y[3] <- mean(sigma_y_IVW_meta)

data1 <- data.frame(bias,EM_var,est_var,cover,est_var_y)


times = 100
reps = 100
TwoStage_est = rep(0,times*reps)
cover_TwoStage_est = rep(0,times*reps)
sigma_TwoStage = rep(0,times*reps)
sigma_y_TwoStage = rep(0,times*reps)
IVW_est = rep(0,times*reps)
cover_IVW_est = rep(0,times*reps)
sigma_IVW = rep(0,times*reps)
sigma_y_IVW = rep(0,times*reps)
IVW_meta_est = rep(0,times*reps)
cover_IVW_meta_est = rep(0,times*reps)
sigma_IVW_meta = rep(0,times*reps)
sigma_y_IVW_meta = rep(0,times*reps)
IVWs_est = rep(0,times*reps)
cover_IVWs_est = rep(0,times*reps)
sigma_IVWs = rep(0,times*reps)
setwd("/spin1/users/zhangh24/MR_MA/")
i2 = 2
total = 0
for(i1 in 1:reps){
  load(paste0("./result/simulation/simulation_",i1,"_",i2,".Rdata"))
  temp = length(result[[2]][[1]])
  TwoStage_est[total+(1:temp)] <- result[[2]][[1]]
  cover_TwoStage_est[total+(1:temp)] = result[[2]][[2]]
  sigma_TwoStage[total+(1:temp)] = result[[2]][[3]]
  sigma_y_TwoStage[total+(1:temp)] = result[[2]][[4]]
  IVW_est[total+(1:temp)] = result[[2]][[5]]
  cover_IVW_est[total+(1:temp)] = result[[2]][[6]]
  sigma_IVW[total+(1:temp)] = result[[2]][[7]]
  sigma_y_IVW[total+(1:temp)] = result[[2]][[8]]
  IVW_meta_est[total+(1:temp)] = result[[2]][[9]]
  cover_IVW_meta_est[total+(1:temp)] = result[[2]][[10]]
  sigma_IVW_meta[total+(1:temp)] = result[[2]][[11]]
  sigma_y_IVW_meta[total+(1:temp)] = result[[2]][[12]]
  IVWs_est[total+(1:temp)] = result[[2]][[13]]
  cover_IVWs_est[total+(1:temp)] = result[[2]][[14]]
  sigma_IVWs[total+(1:temp)] = result[[2]][[15]]
  total = total+temp
}
# IVWs_est_new = IVWs_est[order(IVWs_est)][(times*reps*0.05):(times*reps*0.95)]

bias <- rep(NA,4)
bias[1] <- mean(TwoStage_est)-beta_M
bias[2] <- mean(IVW_est)-beta_M
bias[3] <- mean(IVW_meta_est)-beta_M
bias[4] <- mean(IVWs_est)-beta_M
EM_var <- rep(0,4)
EM_var[1] <- var(TwoStage_est)
EM_var[2] <- var(IVW_est)
EM_var[3] <- var(IVW_meta_est)
EM_var[4] <- var(IVWs_est)
est_var <- rep(0,4)
est_var[1] <- mean(sigma_TwoStage)
est_var[2] <- mean(sigma_IVW)
est_var[3] <- mean(sigma_IVW_meta)
est_var[4] <- mean(sigma_IVWs)
cover <- rep(0,4)
cover[1] <- mean(cover_TwoStage_est)
cover[2] <- mean(cover_IVW_est)
cover[3] <- mean(cover_IVW_meta_est,na.rm=T)
cover[4] <- mean(cover_IVWs_est,na.rm = T)

est_var_y <- rep(NA,4)
est_var_y[1] <- mean(sigma_y_TwoStage)
est_var_y[2] <- mean(sigma_y_IVW)
est_var_y[3] <- mean(sigma_y_IVW_meta)

data2 <- data.frame(bias,EM_var,est_var,cover,est_var_y)

data <- rbind(data1,data2)

write.csv(data,file = "./result/simulation/beta_g_0.01_result.csv")

#goal: merge the simulations datasets
times = 100
reps = 100
TwoStage_est = rep(0,times*reps)
cover_TwoStage_est = rep(0,times*reps)
sigma_TwoStage = rep(0,times*reps)
sigma_y_TwoStage = rep(0,times*reps)
IVW_est = rep(0,times*reps)
cover_IVW_est = rep(0,times*reps)
sigma_IVW = rep(0,times*reps)
sigma_y_IVW = rep(0,times*reps)
IVW_meta_est = rep(0,times*reps)
cover_IVW_meta_est = rep(0,times*reps)
sigma_IVW_meta = rep(0,times*reps)
sigma_y_IVW_meta = rep(0,times*reps)
IVWs_est = rep(0,times*reps)
cover_IVWs_est = rep(0,times*reps)
sigma_IVWs = rep(0,times*reps)
setwd("/spin1/users/zhangh24/MR_MA/")
i2 = 1
total = 0
for(i1 in 1:reps){
  load(paste0("./result/simulation/simulation_",i1,"_",i2,".Rdata"))
  temp = length(result[[1]][[1]])
  TwoStage_est[total+(1:temp)] <- result[[1]][[1]]
  cover_TwoStage_est[total+(1:temp)] = result[[1]][[2]]
  sigma_TwoStage[total+(1:temp)] = result[[1]][[3]]
  sigma_y_TwoStage[total+(1:temp)] = result[[1]][[4]]
  IVW_est[total+(1:temp)] = result[[1]][[5]]
  cover_IVW_est[total+(1:temp)] = result[[1]][[6]]
  sigma_IVW[total+(1:temp)] = result[[1]][[7]]
  sigma_y_IVW[total+(1:temp)] = result[[1]][[8]]
  IVW_meta_est[total+(1:temp)] = result[[1]][[9]]
  cover_IVW_meta_est[total+(1:temp)] = result[[1]][[10]]
  sigma_IVW_meta[total+(1:temp)] = result[[1]][[11]]
  sigma_y_IVW_meta[total+(1:temp)] = result[[1]][[12]]
  IVWs_est[total+(1:temp)] = result[[1]][[13]]
  cover_IVWs_est[total+(1:temp)] = result[[1]][[14]]
  sigma_IVWs[total+(1:temp)] = result[[1]][[15]]
  total = total+temp
}
# IVWs_est_new = IVWs_est[order(IVWs_est)][(times*reps*0.05):(times*reps*0.95)]
beta_M = 0.1

bias <- rep(NA,4)
bias[1] <- mean(TwoStage_est)-beta_M
bias[2] <- mean(IVW_est)-beta_M
bias[3] <- mean(IVW_meta_est)-beta_M
bias[4] <- mean(IVWs_est)-beta_M
EM_var <- rep(0,4)
EM_var[1] <- var(TwoStage_est)
EM_var[2] <- var(IVW_est)
EM_var[3] <- var(IVW_meta_est)
EM_var[4] <- var(IVWs_est)
est_var <- rep(0,4)
est_var[1] <- mean(sigma_TwoStage)
est_var[2] <- mean(sigma_IVW)
est_var[3] <- mean(sigma_IVW_meta)
est_var[4] <- mean(sigma_IVWs)
cover <- rep(0,4)
cover[1] <- mean(cover_TwoStage_est)
cover[2] <- mean(cover_IVW_est)
cover[3] <- mean(cover_IVW_meta_est,na.rm=T)
cover[4] <- mean(cover_IVWs_est,na.rm = T)

est_var_y <- rep(NA,4)
est_var_y[1] <- mean(sigma_y_TwoStage)
est_var_y[2] <- mean(sigma_y_IVW)
est_var_y[3] <- mean(sigma_y_IVW_meta)

data1 <- data.frame(bias,EM_var,est_var,cover,est_var_y)










times = 100
reps = 100
TwoStage_est = rep(0,times*reps)
cover_TwoStage_est = rep(0,times*reps)
sigma_TwoStage = rep(0,times*reps)
sigma_y_TwoStage = rep(0,times*reps)
IVW_est = rep(0,times*reps)
cover_IVW_est = rep(0,times*reps)
sigma_IVW = rep(0,times*reps)
sigma_y_IVW = rep(0,times*reps)
IVW_meta_est = rep(0,times*reps)
cover_IVW_meta_est = rep(0,times*reps)
sigma_IVW_meta = rep(0,times*reps)
sigma_y_IVW_meta = rep(0,times*reps)
IVWs_est = rep(0,times*reps)
cover_IVWs_est = rep(0,times*reps)
sigma_IVWs = rep(0,times*reps)
setwd("/spin1/users/zhangh24/MR_MA/")
i2 = 2
total = 0
for(i1 in 1:reps){
  load(paste0("./result/simulation/simulation_",i1,"_",i2,".Rdata"))
  temp = length(result[[2]][[1]])
  TwoStage_est[total+(1:temp)] <- result[[2]][[1]]
  cover_TwoStage_est[total+(1:temp)] = result[[2]][[2]]
  sigma_TwoStage[total+(1:temp)] = result[[2]][[3]]
  sigma_y_TwoStage[total+(1:temp)] = result[[2]][[4]]
  IVW_est[total+(1:temp)] = result[[2]][[5]]
  cover_IVW_est[total+(1:temp)] = result[[2]][[6]]
  sigma_IVW[total+(1:temp)] = result[[2]][[7]]
  sigma_y_IVW[total+(1:temp)] = result[[2]][[8]]
  IVW_meta_est[total+(1:temp)] = result[[2]][[9]]
  cover_IVW_meta_est[total+(1:temp)] = result[[2]][[10]]
  sigma_IVW_meta[total+(1:temp)] = result[[2]][[11]]
  sigma_y_IVW_meta[total+(1:temp)] = result[[2]][[12]]
  IVWs_est[total+(1:temp)] = result[[2]][[13]]
  cover_IVWs_est[total+(1:temp)] = result[[2]][[14]]
  sigma_IVWs[total+(1:temp)] = result[[2]][[15]]
  total = total+temp
}
# IVWs_est_new = IVWs_est[order(IVWs_est)][(times*reps*0.05):(times*reps*0.95)]

bias <- rep(NA,4)
bias[1] <- mean(TwoStage_est)-beta_M
bias[2] <- mean(IVW_est)-beta_M
bias[3] <- mean(IVW_meta_est)-beta_M
bias[4] <- mean(IVWs_est)-beta_M
EM_var <- rep(0,4)
EM_var[1] <- var(TwoStage_est)
EM_var[2] <- var(IVW_est)
EM_var[3] <- var(IVW_meta_est)
EM_var[4] <- var(IVWs_est)
est_var <- rep(0,4)
est_var[1] <- mean(sigma_TwoStage)
est_var[2] <- mean(sigma_IVW)
est_var[3] <- mean(sigma_IVW_meta)
est_var[4] <- mean(sigma_IVWs)
cover <- rep(0,4)
cover[1] <- mean(cover_TwoStage_est)
cover[2] <- mean(cover_IVW_est)
cover[3] <- mean(cover_IVW_meta_est,na.rm=T)
cover[4] <- mean(cover_IVWs_est,na.rm = T)

est_var_y <- rep(NA,4)
est_var_y[1] <- mean(sigma_y_TwoStage)
est_var_y[2] <- mean(sigma_y_IVW)
est_var_y[3] <- mean(sigma_y_IVW_meta)

data2 <- data.frame(bias,EM_var,est_var,cover,est_var_y)

data <- rbind(data1,data2)

write.csv(data,file = "./result/simulation/beta_g_0.01_result.csv")













result2 = list(TwoStage_est,IVW_est,IVWs_est,
               IVW_est1,
               cover_TwoStage_est,cover_IVW_est,
               cover_IVWs_est,
               cover_IVW_est1,
               sigma_TwoStage,
               sigma_IVW,
               sigma_IVWs,
               sigma_IVW1,
               twostage.nsnps,
               twostage.prop,
               TwoStage_est_all,
               IVW_est_all,
               IVWs_est_all,
               IVW_est_all1)
total = 0
n_thres = 7
TwoStage_est_all <- matrix(0,times*reps,n_thres)
IVW_est_all <- matrix(0,times*reps,n_thres)
IVWs_est_all <- matrix(0,times*reps,n_thres)
IVW_est_all1 <- matrix(0,times*reps,n_thres)
total = 0
for(i1 in 1:reps){
  load(paste0("./result/simulation/simulation_",i1,".Rdata"))
  temp = length(result[[3]][[1]])
  TwoStage_est[total+(1:temp)] <- result[[2]][[1]]
  IVW_est[total+(1:temp)] <- result[[2]][[2]]
  IVWs_est[total+(1:temp)] <- result[[2]][[3]]
  IVW_est1[total+(1:temp)] <- result[[2]][[4]]
  cover_TwoStage_est[total+(1:temp)] <- result[[2]][[5]]
  cover_IVW_est[total+(1:temp)] <- result[[2]][[6]]
  cover_IVWs_est[total+(1:temp)] = result[[2]][[7]]
  cover_IVW_est1[total+(1:temp)] <- result[[2]][[8]]
  sigma_TwoStage[total+(1:temp)] = result[[2]][[9]]
  sigma_IVW[total+(1:temp)] = result[[2]][[10]]
  sigma_IVWs[total+(1:temp)] = result[[2]][[11]]
  sigma_IVW1[total+(1:temp)] = result[[2]][[12]]
  twostage.nsnps[total+(1:temp)] = result[[2]][[13]]
  twostage.prop[total+(1:temp)] = result[[2]][[14]]
  TwoStage_est_all[total+(1:temp),] = result[[2]][[15]]
  IVW_est_all[total+(1:temp),] = result[[2]][[16]]
  IVWs_est_all[total+(1:temp),] = result[[2]][[17]]
  IVW_est_all1[total+(1:temp),] = result[[2]][[18]]
  total = total+temp
}
beta_M = 0.1
mean(TwoStage_est)-beta_M
mean(IVW_est)-beta_M
mean(IVWs_est)-beta_M
mean(IVW_est1)-beta_M
mean(cover_TwoStage_est)
mean(cover_IVW_est)
mean(cover_IVWs_est)
mean(cover_IVW_est1)
var(TwoStage_est)
var(IVW_est)
var(IVWs_est)
var(IVW_est1)
mean(sigma_TwoStage)
mean(sigma_IVW)
mean(sigma_IVWs)
mean(sigma_IVW1)
mean(sigma_y_TwoStage)
mean(sigma_y_IVW)
mean(sigma_y_IVW1)
mean(twostage.nsnps,na.rm = T)
mean(twostage.prop,na.rm = T)

TwoStage_est_var = apply(TwoStage_est_all,2,function(x){var(x,na.rm=T)})
IVW_est_all_var = apply(IVW_est_all,2,function(x){var(x,na.rm=T)})
IVW_est_all_var1 = apply(IVW_est_all1,2,function(x){var(x,na.rm=T)})

p_thres = c(5E-04,1E-03,5E-03,1E-02,5E-02,1E-01,5E-01)

data <- data.frame(-log10(p_thres),TwoStage_est_var,
                   IVW_est_all_var,IVW_est_all_var1)
colnames(data) <- c("-log10(p)","Two_stage","IVW_meta","IVW")
library(reshape2)
data_plot = melt(data,id.var = "-log10(p)")
colnames(data_plot) = c("p","method","emprical_var")

library(ggplot2)
ggplot(data_plot)+
  geom_line(aes(x=p,y=emprical_var,color = method))+
  xlab("-log10(p)")

#apply(IVWs_est_all,2,function(x){var(x,na.rm=T)})





total = 0
n_thres = 7
TwoStage_est_all <- matrix(0,times*reps,n_thres)
IVW_est_all <- matrix(0,times*reps,n_thres)
IVWs_est_all <- matrix(0,times*reps,n_thres)
IVW_est_all1 <- matrix(0,times*reps,n_thres)
total = 0
for(i1 in 1:reps){
  load(paste0("./result/simulation/simulation_",i1,".Rdata"))
  temp = length(result[[4]][[1]])
  TwoStage_est[total+(1:temp)] <- result[[4]][[1]]
  IVW_est[total+(1:temp)] <- result[[4]][[2]]
  IVWs_est[total+(1:temp)] <- result[[4]][[3]]
  IVW_est1[total+(1:temp)] <- result[[4]][[4]]
  cover_TwoStage_est[total+(1:temp)] <- result[[4]][[5]]
  cover_IVW_est[total+(1:temp)] <- result[[4]][[6]]
  cover_IVWs_est[total+(1:temp)] = result[[4]][[7]]
  cover_IVW_est1[total+(1:temp)] <- result[[4]][[8]]
  sigma_TwoStage[total+(1:temp)] = result[[4]][[9]]
  sigma_IVW[total+(1:temp)] = result[[4]][[10]]
  sigma_IVWs[total+(1:temp)] = result[[4]][[11]]
  sigma_IVW1[total+(1:temp)] = result[[4]][[12]]
  twostage.nsnps[total+(1:temp)] = result[[4]][[13]]
  twostage.prop[total+(1:temp)] = result[[4]][[14]]
  TwoStage_est_all[total+(1:temp),] = result[[4]][[15]]
  IVW_est_all[total+(1:temp),] = result[[4]][[16]]
  IVWs_est_all[total+(1:temp),] = result[[4]][[17]]
  IVW_est_all1[total+(1:temp),] = result[[4]][[18]]
  total = total+temp
}
beta_M = 0.1
mean(TwoStage_est)-beta_M
mean(IVW_est)-beta_M
mean(IVWs_est)-beta_M
mean(IVW_est1)-beta_M
mean(cover_TwoStage_est)
mean(cover_IVW_est)
mean(cover_IVWs_est)
mean(cover_IVW_est1)
var(TwoStage_est)
var(IVW_est)
var(IVWs_est)
var(IVW_est1)
mean(sigma_TwoStage)
mean(sigma_IVW)
mean(sigma_IVWs)
mean(sigma_IVW1)
mean(sigma_y_TwoStage)
mean(sigma_y_IVW)
mean(sigma_y_IVW1)
mean(twostage.nsnps,na.rm = T)
mean(twostage.prop,na.rm = T)

TwoStage_est_var = apply(TwoStage_est_all,2,function(x){var(x,na.rm=T)})
IVW_est_all_var = apply(IVW_est_all,2,function(x){var(x,na.rm=T)})
IVW_est_all_var1 = apply(IVW_est_all1,2,function(x){var(x,na.rm=T)})

p_thres = c(5E-04,1E-03,5E-03,1E-02,5E-02,1E-01,5E-01)

data <- data.frame(-log10(p_thres),TwoStage_est_var,
                   IVW_est_all_var,IVW_est_all_var1)
colnames(data) <- c("-log10(p)","Two_stage","IVW_meta","IVW")
library(reshape2)
data_plot = melt(data,id.var = "-log10(p)")
colnames(data_plot) = c("p","method","emprical_var")

library(ggplot2)
ggplot(data_plot)+
  geom_line(aes(x=p,y=emprical_var,color = method))+
  xlab("-log10(p)")
