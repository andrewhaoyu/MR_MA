#plot the ratio estimate distribution
n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.00,0.01,0.03,0.05)
beta_vec <- c(0,0.3,0.5,1)
setwd("/Users/zhangh24/GoogleDrive/MR_MA")
times = 100000


#i1 correponding to n
#i2 corresponding to alpha
library(ggplot2)
n.row <- length(alpha_vec)
n.col <- length(n_vec)
ratio_est_list <- list()
ratio_cover_list <- list()
ci_low_ratio_list <- list()
ci_high_ratio_list <- list()
ratio_est_c_list <- list()
ratio_var_c_list <- list()
ratio_cover_c_list <- list()
ci_low_ratio_c_list <- list()
ci_high_ratio_c_list <- list()
ratio_est_AR_list <- list()
ratio_est_AR_low_list <- list()
ratio_est_AR_high_list <- list()
cover_AR_list <- list()
ratio_est_MR_list <- list()
ratio_est_MR_low_list <- list()
ratio_est_MR_high_list <- list()
cover_MR_list <- list()

temp <- 1
load("./result/simulation/IVW/IVW_merged.Rdata")
for(i4 in 1:4){
  ratio_est <- matrix(0,n.row,n.col)
  ratio_cover <- matrix(0,n.row,n.col)
  ci_low_ratio <- matrix(0,n.row,n.col)
  ci_high_ratio <- matrix(0,n.row,n.col)
  ratio_est_c <- matrix(0,n.row,n.col)
  ratio_cover_c <- matrix(0,n.row,n.col)
  ci_low_ratio_c <- matrix(0,n.row,n.col)
  ci_high_ratio_c <- matrix(0,n.row,n.col)
  ci_high_ratio <- matrix(0,n.row,n.col)
  ratio_est_AR <- matrix(0,n.row,n.col)
  ratio_est_AR_low <- matrix(0,n.row,n.col)
  ratio_est_AR_high <- matrix(0,n.row,n.col)
  cover_AR <- matrix(0,n.row,n.col)
  ratio_est_MR <- matrix(0,n.row,n.col)
  ratio_est_MR_low <- matrix(0,n.row,n.col)
  ratio_est_MR_high <- matrix(0,n.row,n.col)
  cover_MR <- matrix(0,n.row,n.col)
  for(i1 in 1:3){
    for(i2 in 1:4){
      # 
      temp = 12*(i4-1)+4*(i1-1)+i2
      n <- n_vec[i1]
      alpha_G = alpha_vec[i2]
      beta_M = beta_vec[i4]
      result <- result_final[[temp]]
      ratio_est[i2,i1] <- mean(result[[5]])
      ci_low_ratio[i2,i1] <- mean(result[[8]])
      ci_high_ratio[i2,i1] <- mean(result[[9]])
      ratio_cover[i2,i1] <- mean(result[[7]])
      ratio_est_c[i2,i1] <- mean(result[[10]])
      ci_low_ratio_c[i2,i1] <- mean(result[[14]])
      ci_high_ratio_c[i2,i1] <- mean(result[[15]])
      ratio_cover_c[i2,i1] <- mean(result[[12]])
      ratio_est_AR[i2,i1] <- mean(result[[16]])
      ratio_est_AR_low[i2,i1] <- mean(result[[17]],na.rm=T)
      ratio_est_AR_high[i2,i1] <- mean(result[[18]],na.rm=T)
      cover_AR[i2,i1] <- mean(result[[19]])
      ratio_est_MR[i2,i1] <- mean(result[[20]])
      ratio_est_MR_low[i2,i1] <- mean(result[[21]])
      ratio_est_MR_high[i2,i1] <- mean(result[[22]])
      cover_MR[i2,i1] <- mean(result[[23]])
      temp <- temp+1
    }
  }
  ratio_est_list[[i4]] <- ratio_est
  ratio_cover_list[[i4]] <- ratio_cover
  ci_low_ratio_list[[i4]] <- ci_low_ratio
  ci_high_ratio_list[[i4]] <- ci_high_ratio
  ratio_cover_c_list[[i4]] <-   ratio_cover_c
  ci_low_ratio_c_list[[i4]] <- ci_low_ratio_c
  ci_high_ratio_c_list[[i4]] <-  ci_high_ratio_c
  ratio_est_AR_list[[i4]] <- ratio_est_AR
  ratio_est_AR_low_list[[i4]] <- ratio_est_AR_low
  ratio_est_AR_high_list[[i4]] <- ratio_est_AR_high
  cover_AR_list[[i4]] <-   cover_AR
  ratio_est_MR_list[[i4]] <- ratio_est_MR
  ratio_est_MR_low_list[[i4]] <- ratio_est_MR_low
  ratio_est_MR_high_list[[i4]] <- ratio_est_MR_high
  cover_MR_list[[i4]] <-   cover_MR
}

ratio_cover_table <- round(rbind(ratio_cover_list[[1]],
                                 ratio_cover_list[[2]],
                                 ratio_cover_list[[3]],
                                 ratio_cover_list[[4]]),2)
write.csv(ratio_cover_table,file = "./result/simulation/IVW/cover_cover_table.csv")
ratio_cover_c_table <- round(rbind(ratio_cover_c_list[[1]],
                                 ratio_cover_c_list[[2]],
                                 ratio_cover_c_list[[3]],
                                 ratio_cover_c_list[[4]]),2)
write.csv(ratio_cover_c_table,file = "./result/simulation/IVW/ratio_cover_c_table.csv")
cover_AR_table <- round(rbind(cover_AR_list[[1]],
                              cover_AR_list[[2]],
                              cover_AR_list[[3]],
                              cover_AR_list[[4]]),2)
write.csv(cover_AR_table,file = "./result/simulation/IVW/cover_AR_table.csv")

cover_MR_table <- round(rbind(cover_MR_list[[1]],
                              cover_MR_list[[2]],
                              cover_MR_list[[3]],
                              cover_MR_list[[4]]),2)
write.csv(cover_MR_table,file = "./result/simulation/IVW/cover_MR_table.csv")

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

