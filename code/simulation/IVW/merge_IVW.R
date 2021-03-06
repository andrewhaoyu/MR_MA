#Merge the ratio estiamtes results
setwd("/data/zhangh24/MR_MA/")
#setwd("/n/holystore01/LABS/xlin/Lab/hzhang/MR_MA")
n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.0,0.01,0.03,0.05)
beta_vec <-  c(0,0.3,0.5,1)

times = 1000*100
replicates <- 100
p <- 5
result_final <- list()
temp.idx <- 1
for(i4 in 1:4){
  for(i1 in 1:3){
    for(i2 in 1:4){
      Gamma_est <- matrix(0,times,p)
      Gamma_var <- matrix(0,times,p)
      gamma_est <- matrix(0,times,p)
      gamma_var <- matrix(0,times,p)
      ratio_est <- rep(0,times)
      ratio_var <- rep(0,times)
      ratio_cover <- rep(0,times)
      ci_low_ratio<- rep(0,times)
      ci_high_ratio <- rep(0,times)
      ratio_est_c <- rep(0,times)
      ratio_var_c <- rep(0,times)
      ratio_cover_c <- rep(0,times)
      ratio_cover_c <- rep(0,times)
      ci_low_ratio_c <- rep(0,times)
      ci_high_ratio_c <- rep(0,times)
      ratio_est_AR <- rep(0,times)
      ratio_AR_low <- rep(0,times)
      ratio_AR_high <- rep(0,times)
      cover_AR <- rep(0,times)
      ratio_est_MR <- rep(0,times)
      ratio_MR_low <- rep(0,times)
      ratio_MR_high <- rep(0,times)
      cover_MR<- rep(0,times)
    
      total <- 0
      for(i3 in 1:100){
        load(paste0("./result/simulation/IVW/ratio_estimate_",i1,"_",i2,"_",i3,"_",i4,".Rdata"))
        temp <- nrow(result[[1]])
        Gamma_est[total+(1:temp),] <- result[[1]]
        Gamma_var[total+(1:temp),] <- result[[2]]
        gamma_est [total+(1:temp),] <- result[[3]]
        gamma_var[total+(1:temp),] <- result[[4]]
        ratio_est[total+(1:temp)] <- result[[5]]
        ratio_var[total+(1:temp)] <- result[[6]]
        ratio_cover[total+(1:temp)] <- result[[7]]
        ci_low_ratio[total+(1:temp)] <- result[[8]]
        ci_high_ratio[total+(1:temp)] <- result[[9]]
        ratio_est_c[total+(1:temp)] <- result[[10]]
        ratio_var_c[total+(1:temp)] <- result[[11]]
        ratio_cover_c[total+(1:temp)] <- result[[12]]
        ci_low_ratio_c[total+(1:temp)] <- result[[14]]
        ci_high_ratio_c[total+(1:temp)] <- result[[15]]
        ratio_est_AR[total+(1:temp)] <- result[[16]]
        ratio_AR_low[total+(1:temp)] <- result[[17]]
        ratio_AR_high[total+(1:temp)] <- result[[18]]
        cover_AR[total+(1:temp)] <- result[[19]]
        ratio_est_MR[total+(1:temp)] <- result[[20]]
        ratio_MR_low[total+(1:temp)] <- result[[21]]
        ratio_MR_high[total+(1:temp)] <- result[[22]]
        cover_MR[total+(1:temp)] <- result[[23]]
        total <- total+temp
        
      }
      #two loops
      #first loop sample size
      #inner loop for alpha_G
      result_final[[temp.idx]] <-  list(
                                        Gamma_est,
                                        Gamma_var,
                                        gamma_est,
                                        gamma_var,
                                        ratio_est,
                                        ratio_var,
                                        ratio_cover,
                                        ci_low_ratio,
                                        ci_high_ratio,
                                        ratio_est_c,
                                        ratio_var_c,
                                        ratio_cover_c,
                                        ratio_cover_c,
                                        ci_low_ratio_c,
                                        ci_high_ratio_c,
                                        ratio_est_AR,
                                        ratio_AR_low,
                                        ratio_AR_high,
                                        cover_AR,
                                        ratio_est_MR,
                                        ratio_MR_low,
                                        ratio_MR_high,
                                        cover_MR
      )
      temp.idx <- temp.idx+1
    }
  }
  
}


save(result_final,file = paste0("./result/simulation/IVW/IVW_merged.Rdata"))


temp <- 1
result <- result_final[[temp]]
Gamma_est = result[[1]]
Gamma_var = result[[2]]
gamma_est = result[[3]]
gamma_var = result[[4]]
ratio_est = result[[5]]
ratio_var = result[[6]]
ratio_cover = result[[7]]
cover_ratio = result[[8]]
cover_true = result[[9]]
cover_epi = result[[10]]
cover_exact = result[[11]]
cover_true_exact = result[[12]]
mean(cover_ratio)
mean(cover_true)
mean(cover_epi)
mean(cover_exact)
mean(cover_true_exact)
idx.temp <- which(cover_true!=cover_epi)

Ratio = function(Gamma,var_Gamma,gamma,var_gamma,n){
  ratio_est = Gamma/gamma
  var_ratio = var_Gamma/gamma^2+var_gamma*Gamma^2/gamma^4
  n.simu <- 1000000
  z_Gamma <- rnorm(n.simu)
  z_gamma <- rnorm(n.simu,mean = gamma*sqrt(n),sd = 1)
  true_distribution <- z_Gamma/sqrt(1+z_Gamma^2/z_gamma^2)
  q_result <- quantile(true_distribution,c(0.025,0.975))
  z_est <- ratio_est/sqrt(var_ratio)
  cover = ifelse(z_est>=q_result[1]&
                   z_est<=q_result[2],1,0)
  ci_low <- ratio_est-q_result[2]*sqrt(var_ratio)
  ci_high <- ratio_est-q_result[1]*sqrt(var_ratio)
  return(c(ratio_est,var_ratio,cover,ci_low,ci_high))  
}
k <- 1
n <- 15000
alpha_G = 0.00
Gamma = Gamma_est[idx.temp[k]]
var_Gamma = Gamma_var[idx.temp[k]]
gamma = gamma_est[idx.temp[k]]
var_gamma = gamma_var[idx.temp[k]]







ci_low_ratio = result[[13]]
ci_high_ratio = result[[14]]
ci_low_epi = result[[15]]
ci_high_epi = result[[16]]
ci_low_exact = result[[17]]
ci_high_exact = result[[18]]




# library(ggplot2)
# data <- data.frame(z_est,standard_norm,true_distribution)
# colnames(data) <- c("Ratio","Standard Normal","Derived_distribution")
# library(reshape2)
# data.m <- melt(data)
# ggplot(data.m,aes(value,colour=variable))+
#   geom_density()+
#   theme_Publication()




