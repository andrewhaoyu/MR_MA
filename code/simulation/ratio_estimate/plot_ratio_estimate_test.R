#plot the ratio estimate distribution
n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.00,0.01,0.03,0.05)
setwd("/Users/zhangh24/GoogleDrive/MR_MA")
times = 100000


#i1 correponding to n
#i2 corresponding to alpha
library(ggplot2)
n.row <- length(alpha_vec)
n.col <- length(n_vec)
ratio_p1 <- matrix(0,n.row,n.col)
ratio_p2 <- matrix(0,n.row,n.col)
p_weight <- matrix(0,n.row,n.col)
p_delta <- matrix(0,n.row,n.col)
p <- list()
p_ratio <- list()
temp <- 1
load("./result/simulation/ratio_estimate/ratio_estimate_merged_test.Rdata")
for(i1 in 1:3){
  for(i2 in 1:4){
    result <- result_final[[temp]]
    Gamma = result[[1]]
    var_Gamma = result[[2]]
    gamma = result[[3]]
    var_gamma = result[[4]]
    ratio_est = result[[5]]
    var_ratio <- result[[6]]
    
    n <- length(Gamma)
    ratio_p1[i2,i1] <- sum(result[[7]]<=0.05)/n
    ratio_p2[i2,i1] <- sum(result[[8]]<=0.05)/n
    p_weight[i2,i1] <- sum(result[[9]]<=0.05)/n
    p_delta[i2,i1] <- sum(result[[10]]<=0.05)/n
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
grid.arrange(p[[1]],p[[5]],p[[9]],
             p[[2]],p[[6]],p[[10]],
             p[[3]],p[[7]],p[[11]],
             p[[4]],p[[8]],p[[12]],
             ncol=3)
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
temp =1 
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
colnames(data) <- c("Proposed method","IVW","Empirical distribution")
library(reshape2)
data.m <- melt(data)
data.m.temp <- data.m

ggplot(data.m,aes(value,colour=variable))+
  geom_density()+
  theme_Publication()+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(face="bold"))+
  scale_fill_discrete(name = "New Legend Title")
dev.off()

