times <- 1000000
n <- 15000
sigma_y=sigma_m = 1
alpha_G_vec <- seq(0.012,0.048,
                   by = (0.048-0.012)/5)
p <- length(alpha_G_vec)
alpha_G = 0.03
alpha_G_value <- rep(0,times*p)
true_distribution_value <- rep(0,times*p)
total <- 0
for(i in 1:p){
  z_Gamma <- rnorm(times,mean = 0, sd = sqrt(sigma_y/n))
  z_gamma <- rnorm(times,mean = alpha_G_vec[i],sd = sqrt(sigma_m/n))
  true_distribution_temp <- z_Gamma/z_gamma
  q_temp <- quantile(true_distribution_temp,c(0.025,0.975))
  idx1 = which(true_distribution_temp>=q_temp[1]&
                 true_distribution_temp<=q_temp[2])
  true_distribution_value[total+(1:length(idx1))] <- true_distribution_temp[idx1]
  alpha_G_value[total+(1:length(idx1))] <- alpha_G_vec[i]
  total = total+length(idx1)
}
true_distribution_value <- true_distribution_value[1:total]
alpha_G_value <- alpha_G_value[1:total]
library(ggplot2)
new.data <- data.frame(alpha_G_value=paste0("alpha=",alpha_G_value),
                      true_distribution_value)
ggplot(new.data,aes(true_distribution_value,
                    color=alpha_G_value))+
  geom_density()
