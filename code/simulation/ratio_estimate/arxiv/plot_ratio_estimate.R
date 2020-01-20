#plot the ratio estimate distribution
n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.01,0.03,0.05)
setwd("/Users/zhangh24/GoogleDrive/MR_MA")
times = 100000


#i1 correponding to n
#i2 corresponding to alpha
library(ggplot2)
cover_ratio <- matrix(0,3,3)
cover_true <- matrix(0,3,3)
cover_epi <- matrix(0,3,3)
cover_exact <- matrix(0,3,3)
cover_true_exact <- matrix(0,3,3)
ci_low_ratio <- matrix(0,3,3)
ci_high_ratio <- matrix(0,3,3)
ci_ratio <- matrix(0,3,3)
ci_low_epi <- matrix(0,3,3)
ci_high_epi <- matrix(0,3,3)
ci_epi <- matrix(0,3,3)
ci_exact <- matrix(0,3,3)
ci_low_exact <- matrix(0,3,3)
ci_high_exact <- matrix(0,3,3)

p <- list()
p_ratio <- list()
temp <- 1
load("./result/simulation/ratio_estimate/ratio_estimate_merged.Rdata")
for(i1 in 1:3){
  for(i2 in 1:3){
    result <- result_final[[temp]]
    Gamma = result[[1]]
    var_Gamma = result[[2]]
    gamma = result[[3]]
    var_gamma = result[[4]]
    var_ratio <- result[[6]]
    n <- length(Gamma)
    cover_ratio[i2,i1] <- mean(result[[8]])
    cover_true[i2,i1] <- mean(result[[9]])
    cover_epi[i2,i1] <- mean(result[[10]])
    cover_exact[i2,i1] <- mean(result[[11]])
    cover_true_exact[i2,i1] <- mean(result[[12]])
    ci_low_ratio[i2,i1] <- mean(result[[13]])
    ci_high_ratio[i2,i1] <- mean(result[[14]])
    ci_ratio[i2,i1] <- paste0(ci_low_ratio[i2,i1],", ",ci_high_ratio[i2,i1])
    ci_low_epi[i2,i1] <- mean(result[[15]])
    ci_high_epi[i2,i1] <- mean(result[[16]])
    ci_epi[i2,i1] <- paste0(ci_low_epi[i2,i1],", ",ci_high_epi[i2,i1])
    ci_low_exact[i2,i1] <- mean(result[[17]])
    ci_high_exact[i2,i1] <- mean(result[[18]])
    ci_exact[i2,i1] <- paste0(ci_low_exact[i2,i1],", ",ci_high_exact[i2,i1])
ratio_est = result[[5]]
ratio_var = result[[6]]
z_est = ratio_est/sqrt(ratio_var)
standard_norm = rnorm(times)
z_Gamma <- rnorm(times)
z_gamma <- rnorm(times,mean = alpha_vec[i2]*sqrt(n_vec[i1]),sd = 1)

true_distribution <- z_Gamma/sqrt(1+z_Gamma^2/z_gamma^2)

data <- data.frame(z_est,standard_norm,true_distribution)
colnames(data) <- c("Ratio/sd","Standard Normal","Derived_distribution")
library(reshape2)
data.m <- melt(data)
data.m.temp <- data.m
p[[temp]] <- ggplot(data.m,aes(value,colour=variable))+
  geom_density()+
  theme_Publication()+
  theme(legend.position = "none")
quantemp <- quantile(ratio_est,c(0.025,0.975))
idx <- which(ratio_est>=quantemp[1]&
               ratio_est<=quantemp[2])
ratio_new <- ratio_est[idx]
standard_norm <- rnorm(length(idx))


data <- data.frame(ratio_new,standard_norm)
colnames(data) <- c("Ratio","Standard Normal ")
data.m <- melt(data)

p_ratio[[temp]] <- ggplot(data.m,aes(value,colour=variable))+
  geom_density()+
  theme_Publication()+
  theme(legend.position = "none")
temp <- temp+1
  }
}
write.csv(cover_ratio,file = "./result/simulation/ratio_estimate/cover_ratio.csv")
write.csv(cover_true,file = "./result/simulation/ratio_estimate/cover_true.csv")
write.csv(cover_epi,file = "./result/simulation/ratio_estimate/cover_epi.csv")
write.csv(cover_true_exact,file = "./result/simulation/ratio_estimate/cover_true_exact.csv")
write.csv(cover_exact,file = "./result/simulation/ratio_estimate/cover_exact.csv")
write.csv(ci_ratio,file = "./result/simulation/ratio_estimate/ci_ratio.csv")
write.csv(ci_epi,file = "./result/simulation/ratio_estimate/ci_epi.csv")
write.csv(ci_exact,file = "./result/simulation/ratio_estimate/ci_exact.csv")


library(gridExtra)
png("./result/simulation/ratio_estimate/ratio_sd_plot.png",width = 16,height = 8,
    unit = "in",res = 300)
grid.arrange(p[[1]],p[[4]],p[[7]],
             p[[2]],p[[5]],p[[8]],
             p[[3]],p[[6]],p[[9]],ncol=3)
dev.off()

png("./result/simulation/ratio_estimate/ratio_sd_plot_legend.png",width = 8,height = 8,
    unit = "in",res = 300)
ggplot(data.m.temp,aes(value,colour=variable))+
  geom_density()+
  theme_Publication()
dev.off()

png("./result/simulation/ratio_estimate/ratio_plot.png",width = 16,height = 8,
    unit = "in",res = 300)
grid.arrange(p_ratio[[1]],p_ratio[[4]],p_ratio[[7]],
             p_ratio[[2]],p_ratio[[5]],p_ratio[[8]],
             p_ratio[[3]],p_ratio[[6]],p_ratio[[9]],ncol=3)
dev.off()

png("./result/simulation/ratio_estimate/ratio_plot_legend.png",width = 8,height = 8,
    unit = "in",res = 300)
ggplot(data.m,aes(value,colour=variable))+
  geom_density()+
  theme_Publication()
dev.off()
