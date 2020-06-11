setwd("/Users/zhangh24/GoogleDrive/MR_MA/result/simulation/PRS")
load("summary_gwas_06.rdata")
alpha_est <- result_new[[1]]  
gamma_est <- result_new[[2]]
l <- 1
plot(alpha_est[,l],gamma_est[,l])
abline(a=0,b=0.15)
points(0,0,col="red")
model <- lm(gamma_est[,l]~alpha_est[,l])
summary(model)

summary(model)