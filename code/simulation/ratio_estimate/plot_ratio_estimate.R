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
cover_ratio_list <- list()
cover_true_list <- list()
cover_epi_list <- list()
cover_exact_list <- list()
cover_true_exact_list <- list()
ci_low_ratio_list <- list()
ci_high_ratio_list <- list()
ci_ratio_list <- list()
ci_low_epi_list <- list()
ci_high_epi_list <- list()
ci_epi_list <- list()
ci_exact_list <- list()
ci_low_exact_list <- list()
ci_high_exact_list <- list()

p <- list()
p_ratio <- list()
temp <- 1
load("./result/simulation/ratio_estimate/ratio_estimate_merged.Rdata")
for(i4 in 1:4){
  cover_ratio <- matrix(0,n.row,n.col)
  cover_true <- matrix(0,n.row,n.col)
  cover_epi <- matrix(0,n.row,n.col)
  cover_exact <- matrix(0,n.row,n.col)
  cover_true_exact <- matrix(0,n.row,n.col)
  ci_low_ratio <- matrix(0,n.row,n.col)
  ci_high_ratio <- matrix(0,n.row,n.col)
  ci_ratio <- matrix(0,n.row,n.col)
  ci_low_epi <- matrix(0,n.row,n.col)
  ci_high_epi <- matrix(0,n.row,n.col)
  ci_epi <- matrix(0,n.row,n.col)
  ci_exact <- matrix(0,n.row,n.col)
  ci_low_exact <- matrix(0,n.row,n.col)
  ci_high_exact <- matrix(0,n.row,n.col)
  alpha_U = 0.1
  beta_U <- 0.1
  sigma_y = 1
  sigma_m = 1
  for(i1 in 1:3){
    for(i2 in 1:4){
      # 
      # temp = 12*(i4-1)+4*(i1-1)+i2
      n <- n_vec[i1]
      alpha_G = alpha_vec[i2]
      beta_M = beta_vec[i4]
      result <- result_final[[temp]]
      Gamma = result[[1]]
      var_Gamma = result[[2]]
      gamma = result[[3]]
      var_gamma = result[[4]]
      ratio_est = result[[5]]
      var_ratio <- result[[6]]
      n <- length(Gamma)
      cover_ratio[i2,i1] <- mean(result[[7]])
      cover_true[i2,i1] <- mean(result[[8]])
      cover_epi[i2,i1] <- mean(result[[9]])
      cover_exact[i2,i1] <- mean(result[[10]])
      cover_true_exact[i2,i1] <- mean(result[[11]])
      ci_low_ratio[i2,i1] <- mean(result[[12]])
      ci_high_ratio[i2,i1] <- mean(result[[13]])
      ci_ratio[i2,i1] <- paste0(ci_low_ratio[i2,i1],", ",ci_high_ratio[i2,i1])
      ci_low_epi[i2,i1] <- mean(result[[14]])
      ci_high_epi[i2,i1] <- mean(result[[15]])
      ci_epi[i2,i1] <- paste0(round(ci_low_epi[i2,i1],2),", ",round(ci_high_epi[i2,i1],2))
      ci_low_exact[i2,i1] <- mean(result[[16]])
      ci_high_exact[i2,i1] <- mean(result[[17]])
      ci_exact[i2,i1] <- paste0(round(ci_low_exact[i2,i1],2),", ",round(ci_high_exact[i2,i1],2))
    
      true_distribution = ratio_est
      dot <- seq(0.01,0.99,by=0.002)
      q_true <- quantile(true_distribution,dot)
      pdot <-ifelse(dot<=0.5,2*dot,2*(1-dot))
      
      
      n.simu = 1e6
      z_Gamma <- rnorm(n.simu,mean =alpha_G*beta_M,sd =sqrt((sigma_y+beta_M^2*sigma_m+(beta_M*alpha_U)^2)/n))
      z_gamma <- rnorm(n.simu,mean = alpha_G,sd = sqrt((alpha_U^2+sigma_m)/n))
      true_distribution2 <- z_Gamma/z_gamma
        dot2 = dot
        for(i in 1:length(dot)){
          dot2[i] <- sum(true_distribution2<=q_true[i])/length(true_distribution2)
        
        }
        pdot2 <- ifelse(dot2<=0.5,2*dot2,2*(1-dot2))
        data <- data.frame(pdot[order(pdot)],
                           pdot2[order(pdot2)])
        colnames(data) <- c("pdot",
                            "pdot2")
        ggplot(data,aes(pdot,pdot2))+
          geom_point()+
          geom_abline(slope=1)+
          xlab(paste0("Empirical distribution"))+
          ylab(paste0("Derived distribution"))+
          ggtitle(paste0("n=",n,", beta = ",beta_M,", alpha =",alpha_G))+
          theme_Publication()+
          scale_y_continuous(limits=c(0,1))
        
        n.simu = 1e3
        true_distribution3 = rep(0,n.simu*length(Gamma))
        total <- 0
        for(k in 1:length(Gamma)){
          z_Gamma <- rnorm(n.simu,mean =Gamma[k],sd =sqrt(var_Gamma[k]))
          z_gamma <- rnorm(n.simu,mean = gamma[k],sd = sqrt(var_gamma[k]))
          temp_simulation <- z_Gamma/z_gamma  
          true_distribution3[total+(1:n.simu)] <- temp_simulation
          total = total+n.simu
        }
        
        dot3 = dot
        for(i in 1:length(dot)){
          print(i)
          dot3[i] =  sum(true_distribution3<=q_true[i])/length(true_distribution3)
        }
        pdot3 <- ifelse(dot3<=0.5,2*dot3,2*(1-dot3))
        data <- data.frame(pdot[order(pdot)],
                           pdot3[order(pdot3)])
        colnames(data) <- c("pdot",
                            "pdot2")
        ggplot(data,aes(pdot,pdot2))+
          geom_point()+
          geom_abline(slope=1)+
          xlab(paste0("Empirical distribution"))+
          ylab(paste0("Derived distribution"))+
          ggtitle(paste0("n=",n,", beta = ",beta_M,", alpha =",alpha_G))+
          theme_Publication()+
          scale_y_continuous(limits=c(0,1))
        
        
      #     standard_norm = rnorm(times)
      #     z_Gamma <- rnorm(times)
      #     z_gamma <- rnorm(times,mean = alpha_vec[i2]*sqrt(n_vec[i1]),sd = 1)
      #     
      #     true_distribution <- z_Gamma/sqrt(1+z_Gamma^2/z_gamma^2)
      # data <- data.frame(z_est,standard_norm,true_distribution)
      # colnames(data) <- c("Derived distribution","IVW","Real distribution")
      # library(reshape2)
      # data.m <- melt(data)
      # data.m.temp <- data.m
      # p[[temp]] <- ggplot(data.m,aes(value,colour=variable))+
      #   geom_density()+
      #   theme_Publication()+
      #   theme(legend.position = "none")
      # quantemp <- quantile(ratio_est,c(0.025,0.975))
      # idx <- which(ratio_est>=quantemp[1]&
      #                ratio_est<=quantemp[2])
      # ratio_new <- ratio_est[idx]
      # standard_norm <- rnorm(length(idx))
      # 
      # 
      # data <- data.frame(ratio_new,standard_norm)
      # colnames(data) <- c("Ratio","Standard Normal ")
      # data.m <- melt(data)
      # 
      # p_ratio[[temp]] <- ggplot(data.m,aes(value,colour=variable))+
      #   geom_density()+
      #   theme_Publication()+
      #   theme(legend.position = "none")
      temp <- temp+1
    }
  }
  cover_ratio_list[[i4]] <- cover_ratio
  cover_true_list[[i4]] <-   cover_true
  cover_epi_list[[i4]] <- cover_epi
  cover_exact_list[[i4]] <- cover_exact
  cover_true_exact_list[[i4]] <- cover_true
  ci_low_ratio_list[[i4]] <- ci_low_ratio
  ci_high_ratio_list[[i4]] <- ci_high_ratio
  ci_ratio_list[[i4]] <- ci_ratio
  ci_low_epi_list[[i4]] <- ci_low_epi
  ci_high_epi_list[[i4]] <- ci_high_epi
  ci_epi_list[[i4]] <- ci_epi
  ci_exact_list[[i4]] <- ci_exact
  ci_low_exact_list[[i4]] <- ci_low_exact
  ci_high_exact_list[[i4]] <- ci_high_exact
}

cover_epi_table <- round(rbind(cover_epi_list[[1]],
                               cover_epi_list[[2]],
                               cover_epi_list[[3]],
                               cover_epi_list[[4]]),2)
write.csv(cover_epi_table,file = "./result/simulation/ratio_estimate/cover_epi_table.csv")

cover_exact_table <- round(rbind(cover_exact_list[[1]],
                                 cover_exact_list[[2]],
                                 cover_exact_list[[3]],
                                 cover_exact_list[[4]]),2)
write.csv(cover_exact_table,file = "./result/simulation/ratio_estimate/cover_exact_table.csv")

ci_exact_table <- rbind(ci_exact_list[[1]],
                              ci_exact_list[[2]],
                              ci_exact_list[[3]],
                              ci_exact_list[[4]])
write.csv(ci_exact_table,file = "./result/simulation/ratio_estimate/ci_exact_table.csv")

write.csv(cover_ratio,file = "./result/simulation/ratio_estimate/cover_ratio.csv")
write.csv(cover_true,file = "./result/simulation/ratio_estimate/cover_true.csv")
write.csv(cover_epi,file = "./result/simulation/ratio_estimate/cover_epi.csv")
write.csv(cover_true_exact,file = "./result/simulation/ratio_estimate/cover_true_exact.csv")
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

