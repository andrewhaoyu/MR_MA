#plot the ratio estimate distribution
n_vec <- c(15000,75000,150000)
beta_vec <- c(0.01,0.03,0.05)
setwd("/Users/zhangh24/GoogleDrive/MR_MA")
times = 50000
library(ggplot2)
cover_ratio <- matrix(0,3,3)
cover_true <- matrix(0,3,3)
cover_epi <- matrix(0,3,3)
p <- list()
temp <- 1
for(i1 in 1:3){
  for(i2 in 1:3){
    load(paste0("./result/simulation/ratio_estimate/ratio_estimate_",i1,"_",i2,".Rdata"))
    cover_ratio[i1,i2] <- result[[8]]
    cover_true[i1,i2] <- result[[9]]
    cover_epi[i1,i2] <- result[[10]]
ratio_est = result[[5]]
ratio_var = result[[6]]
z_est = ratio_est/sqrt(ratio_var)
standard_norm = rnorm(times)
z_Gamma <- rnorm(times)
z_gamma <- rnorm(times,mean = beta_vec[i2]*sqrt(n_vec[i1]),sd = 1)

true_distribution <- z_Gamma/sqrt(1+z_Gamma^2/z_gamma^2)

data <- data.frame(z_est,standard_norm,true_distribution)
colnames(data) <- c("Ratio","Standard Normal","Derived_distribution")
library(reshape2)
data.m <- melt(data)
p[[temp]] <- ggplot(data.m,aes(value,colour=variable))+
  geom_density()+
  theme_Publication()+
  theme(legend.position = "none")
temp <- temp+1
  }
}

library(gridExtra)
png("./result/simulation/ratio_estimate/ratio_plot.png",width = 8,height = 8,
    unit = "in",res = 300)
grid.arrange(p[[1]],p[[2]],p[[3]],
             p[[4]],p[[5]],p[[6]],
             p[[7]],p[[8]],p[[9]],nrow=3)
dev.off()

png("./result/simulation/ratio_estimate/ratio_plot_legend.png",width = 8,height = 8,
    unit = "in",res = 300)
ggplot(data.m,aes(value,colour=variable))+
  geom_density()+
  theme_Publication()
dev.off()
