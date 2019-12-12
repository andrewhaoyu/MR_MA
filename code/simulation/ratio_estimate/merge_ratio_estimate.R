#Merge the ratio estiamtes results
setwd("/data/zhangh24/MR_MA/")

n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.01,0.03,0.05)

times = 1000*100
replicates <- 100
result_final <- list()
temp.idx <- 1
for(i1 in 1:3){
  for(i2 in 1:3){
     
    Gamma_est <- rep(0,times)
    Gamma_var <- rep(0,times)
    gamma_est <- rep(0,times)
    gamma_var <- rep(0,times)
    ratio_est <- rep(0,times)
    ratio_var <- rep(0,times)
    ratio_cover <- rep(0,times)
    cover_ratio <- rep(0,replicates)
    cover_true <- rep(0,replicates)
    cover_epi <- rep(0,replicates)
    ci_low_ratio <- rep(0,times)
    ci_high_ratio <- rep(0,times)
    ci_low_epi <- rep(0,times)
    ci_high_epi <- rep(0,times)
    total <- 0
    for(i3 in 1:100){
      load(paste0("./result/simulation/ratio_estimate/ratio_estimate_",i1,"_",i2,"_",i3,".Rdata"))
      temp <- length(result[[1]])
      Gamma_est[total+(1:temp)] <- result[[1]]
      Gamma_var[total+(1:temp)] <- result[[2]]
      gamma_est [total+(1:temp)] <- result[[3]]
      gamma_var[total+(1:temp)] <- result[[4]]
      ratio_est[total+(1:temp)] <- result[[5]]
      ratio_var[total+(1:temp)] <- result[[6]]
      ratio_cover[total+(1:temp)] <- result[[7]]
      #cover_ratio cover_true cover_epi already take the mean
      cover_ratio[i3] <- result[[8]]
      cover_true[i3] <- result[[9]]
      cover_epi[i3] <- result[[10]]
      ci_low_ratio[total+(1:temp)] <- result[[11]]
      ci_high_ratio[total+(1:temp)] <- result[[12]]
      ci_low_epi[total+(1:temp)] <- result[[13]]
      ci_high_epi[total+(1:temp)] <- result[[14]]
      total <- total+temp
    }
    #two loops
    #first loop sample size
    #inner loop for alpha_G
    result_final[[temp.idx]] <- list(Gamma_est,
                   Gamma_var,
                   gamma_est,
                   gamma_var,
                   ratio_est,
                   ratio_var,
                   ratio_cover,
                   cover_ratio,
                   cover_true,
                   cover_epi,
                   ci_low_ratio,
                   ci_high_ratio,
                   ci_low_epi,
                   ci_high_epi)
    temp.idx <- temp.idx+1
  }
}
save(result_final,file = paste0("./result/simulation/ratio_estimate/ratio_estimate_merged.Rdata"))





# library(ggplot2)
# data <- data.frame(z_est,standard_norm,true_distribution)
# colnames(data) <- c("Ratio","Standard Normal","Derived_distribution")
# library(reshape2)
# data.m <- melt(data)
# ggplot(data.m,aes(value,colour=variable))+
#   geom_density()+
#   theme_Publication()




