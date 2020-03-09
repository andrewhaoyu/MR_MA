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

