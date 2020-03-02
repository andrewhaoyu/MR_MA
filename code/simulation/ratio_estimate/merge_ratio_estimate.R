#Merge the ratio estiamtes results
setwd("/data/zhangh24/MR_MA/")

n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.0,0.01,0.03,0.05)

times = 1000*100
replicates <- 100
result_final <- list()
temp.idx <- 1
for(i1 in 1:3){
  for(i2 in 1:4){
     
    Gamma_est <- rep(0,times)
    Gamma_var <- rep(0,times)
    gamma_est <- rep(0,times)
    gamma_var <- rep(0,times)
    ratio_est <- rep(0,times)
    ratio_var <- rep(0,times)
    ratio_cover <- rep(0,times)
    cover_ratio <- rep(0,times)
    cover_true <- rep(0,times)
    cover_epi <- rep(0,times)
    cover_exact <- rep(0,times)
    ci_low_ratio <- rep(0,times)
    ci_high_ratio <- rep(0,times)
    ci_low_epi <- rep(0,times)
    ci_high_epi <- rep(0,times)
    ci_low_exact <- rep(0,times)
    ci_high_exact <- rep(0,times)
    
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
      #this ratio_cover is not useful
      #just put there for now
      ratio_cover[total+(1:temp)] <- result[[7]]
      cover_ratio[total+(1:temp)] <- result[[8]]
      cover_true[total+(1:temp)] <- result[[9]]
      cover_epi[total+(1:temp)] <- result[[10]]
      cover_exact[total+(1:temp)] <- result[[11]]
      cover_true_exact[total+(1:temp)] <- result[[12]]
      ci_low_ratio[total+(1:temp)] <- result[[13]]
      ci_high_ratio[total+(1:temp)] <- result[[14]]
      ci_low_epi[total+(1:temp)] <- result[[15]]
      ci_high_epi[total+(1:temp)] <- result[[16]]
      ci_low_exact[total+(1:temp)] <- result[[17]]
      ci_high_exact[total+(1:temp)] <- result[[18]]
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
                                     cover_exact,
                                     cover_true_exact,
                                     ci_low_ratio,
                                     ci_high_ratio,
                                     ci_low_epi,
                                     ci_high_epi,
                                     ci_low_exact,
                                     ci_high_exact)
    temp.idx <- temp.idx+1
  }
}
save(result_final,file = paste0("./result/simulation/ratio_estimate/ratio_estimate_merged.Rdata"))


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




