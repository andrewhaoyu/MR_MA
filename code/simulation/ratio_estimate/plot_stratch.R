#test plot
n.simu <- 1000000
sigma_y = sigma_m = 1
n_vec <- c(15000,75000,150000)
alpha_vec <- c(0.00,0.01,0.03,0.05)
i1 = 1
i2 = 2
n = n_vec[i1]
alpha_G = alpha_vec[i2]
z_Gamma <- rnorm(n.simu,mean = 0, sd = sqrt(sigma_y/n))
z_gamma <- rnorm(n.simu,mean = alpha_G,sd = sqrt(sigma_m/n))
true_distribution <- z_Gamma/z_gamma
q_temp=  quantile(true_distribution,c(0.025,0.975))


z_Gamma <- rnorm(n.simu,mean = 0, sd = 1)
z_gamma <- rnorm(n.simu,mean = 1.224745,sd = 1)
true_distribution <- z_Gamma/z_gamma
q_temp=  quantile(true_distribution,c(0.01,0.99))
idx <- which(true_distribution>=q_temp[1]&
               true_distribution<=q_temp[2])
true_distribution = true_distribution[idx]


i1 = 1
i2 = 1
n = n_vec[i1]

alpha_G = alpha_vec[i2]
z_Gamma <- rnorm(n.simu,mean = 0, sd = sqrt(sigma_y/n))
z_gamma <- rnorm(n.simu,mean = alpha_G,sd = sqrt(sigma_m/n))
true_distribution2 <- z_Gamma/z_gamma
q_temp=  quantile(true_distribution2,c(0.01,0.99))
idx <- which(true_distribution2>=q_temp[1]&
               true_distribution2<=q_temp[2])
true_distribution2 = true_distribution2[idx]

distribution = c(true_distribution,true_distribution2)
alpha = c(rep("0.01",length(true_distribution)),rep("0.00",length(true_distribution)))
data = data.frame(alpha,distribution)
library(ggplot2)
ggplot(data,aes(x=sqrt(n)*distribution,color=alpha))+
  geom_density()


q_high = ratio_est-ci_low_exact
q_low = ratio_est-ci_high_exact
jdx <- which(q_low<=q_result_final[1]|q_high>=q_result_final[2])
idx <- which(ratio_est>=qcauchy(0.025)&
               ratio_est<=qcauchy(0.975))


library(gridExtra)
library(ggplot2)
n.simu <- 1000000
# z_Gamma <- rnorm(n.simu,mean =0,sd =sqrt(var_Gamma))
# z_gamma <- rnorm(n.simu,mean = gamma,sd = sqrt(var_gamma))
#z_Gamma <- rnorm(n.simu,mean =0,sd =sqrt(var_Gamma))
#z_gamma <- rnorm(n.simu,mean = gamma,sd = sqrt(var_gamma))
#z_gamma <- rnorm(n.simu,mean = sqrt(n)*alpha_G,sd = sqrt((n-1)*var_gamma))
sigma_yg_est = 1
sigma_m_est = 1
gamma_vec <- c(0,0.01,0.03,0.05)
for(j in 1:length(gamma_vec)){
  gamma = gamma_vec[j]
  z_Gamma <- rnorm(n.simu,mean =0,sd =sqrt(sigma_yg_est/n))
  z_gamma <- rnorm(n.simu,mean = gamma,sd = sqrt(sigma_m_est/n))
  true_distribution <- z_Gamma/z_gamma
  dot <- seq(0.01,0.99,by=0.001)
  q_true <- quantile(true_distribution,dot)
  pdot <-ifelse(dot<=0.5,2*dot,2*(1-dot))
  
  
  
  
  
  alpha_high = gamma+1.96*sqrt(1/15000)
  alpha_low = gamma-1.96*sqrt(1/15000)
  
  alpha_est_vec = seq(from = alpha_low,to = alpha_high,by = (alpha_high-alpha_low)/(3))
  p<- list()
  for(k in 1:length(alpha_est_vec)){
    alpha_est = alpha_est_vec[k]
    z_Gamma <- rnorm(n.simu,mean =0,sd =sqrt(sigma_yg_est/n))
    z_gamma <- rnorm(n.simu,mean = alpha_est,sd = sqrt(sigma_m_est/n))
    true_distribution2 =  z_Gamma/z_gamma
    dot2 = dot
    for(i in 1:length(dot)){
      idx <- which(true_distribution2<=q_true[i])
      dot2[i] = length(idx)/length(true_distribution2)
    }
    pdot2 <- ifelse(dot2<=0.5,2*dot2,2*(1-dot2))
    data <- data.frame(pdot[order(pdot)],
                       pdot2[order(pdot2)])
    colnames(data) <- c("pdot",
                        "pdot2")
    
    p[[k]] <- ggplot(data,aes(pdot,pdot2))+
      geom_point()+
      geom_abline(slope=1)+
      xlab(paste0("true p-value (alpha = ",gamma,")"))+
      ylab(paste0("estimated p-value ( alpha_hat = ",round(alpha_est,3),")"))+
      theme_Publication()+
      scale_y_continuous(limits=c(0,1))
    
    
  }
  png(filename = paste0("./result/simulation/ratio_estimate/qq_plot_true_plug/alpha_",gamma,".png"),res=300,width=12,height = 6,unit="in")
  grid.arrange(p[[1]],p[[2]],p[[3]],p[[4]],
               nrow=1)
  dev.off()
}








