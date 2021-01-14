#generate phenotypes for the two-sample MR analysis


load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata")
n.sub <- nrow(genotype_s)
n.snp = 5000
beta_M = 0.15
sigma_G  = 0.4
#sigma_m = 1
sigma_m = 1-sigma_G
#generate M phenotypes
set.seed(666)
alpha_G = rnorm(n.snp,sd=sqrt(sigma_G/n.snp))
#alpha_G[1] =0.1
n.rep = 1000
M_mat <- matrix(0,n.sub,n.rep)
#G_value = genotype_s%*%alpha_G + alpha_U*U
G_value = genotype_s%*%alpha_G 
for(j in 1:n.rep){
  print(j)
  M_mat[,j] = G_value+rnorm(n.sub,sd = sqrt(sigma_m))
  #M_mat[,j] = G_value
  #+rnorm(n.sub,sd = sqrt(sigma_e))
    }
save(M_mat,file = "/data/zhangh24/MR_MA/result/simulation/prs/M_mat.rdata")



#generate Y phenotypes
load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_Y.rdata")
n.sub <- nrow(genotype_s2)
n.rep = 1000
Y_mat <- matrix(0,n.sub,n.rep)
G_value2 =  genotype_s2%*%alpha_G
rho = 0.3
sigma_y = 1-beta_M^2
sigma_ym = sigma_y*sigma_m*rho
Sigma = matrix(c(sigma_y,sigma_ym,sigma_ym,sigma_m),2,2)
library(MASS)

for(j in 1:n.rep){
  print(j)
  error = mvrnorm(n.sub,mu = c(0,0),Sigma = Sigma)
  M =G_value2+error[,2]
  #+ rnorm(n.sub,sd = sqrt(sigma_e))  
  #M =G_value2 
  #+ rnorm(n.sub,sd = sqrt(sigma_e))  
  #Y_mat[,j] = beta_M*M+rnorm(n.sub,sd=sqrt(sigma_ey))  
  Y_mat[,j] = beta_M*M+error[,1]
  #+rnorm(n.sub,sd=sqrt(sigma_ey))  
  #Y_mat[,j] = beta_M*M+beta_U*U2+rnorm(n.sub,sd=sqrt(sigma_ey))  
}
save(Y_mat,file = "/data/zhangh24/MR_MA/result/simulation/prs/Y_mat.rdata")

# n.rep = 1
# n.snp = 1000
# beta_est = rep(0,n.snp)
# alpha_est = rep(0,n.snp)
# Gamma_est = rep(0,n.snp)
# library(RcppArmadillo)
# FitLinearmodel <- function(y,x){
#   model <- fastLm(X=x,y=y)
#   if(is.na(coef(model)[2])){
#     result <- c(0,1,1)
#   }else{
#     result <- coef(summary(model))[2,c(1,2,4)]  
#   }
#   return(result)
# }
# n.train = 100000
# M_mat_train = M_mat[(1:n.train),]
# genotype_m_train = genotype_s[(1:n.train),]
# 
# n.snp = 5000
# for(k in 1:n.snp){
#   print(k)
#   #model1 <- lm(M_mat_train[,2]~genotype_m_train[,k])
#   #alpha_est[k] = coef(summary(model1))[2,1]
#   alpha_est[k] = FitLinearmodel(M_mat_train[,1],cbind(1,genotype_m_train[,k]))[1]
#    
#   #model2 <- lm(Y_mat[,2]~genotype_s2[,k])
#   #Gamma_est[k] = coef(summary(model2))[2,1]
#   Gamma_est[k] = FitLinearmodel(Y_mat[,1],cbind(1,genotype_s2[,k]))[1]
#   #beta_est[k] = Gamma_est/alpha_est
# }
# 
# genotype_m_test = genotype_s[(n.train+1):nrow(genotype_s),1:n.snp]
# 
# 
# prs_m_mat <- genotype_s%*%alpha_est
# prs_y_mat <- genotype_s%*%Gamma_est
# model = lm(prs_y_mat~prs_m_mat)
# summary(model)
# 
# 
# 
# 
# 
# prs_y = genotype_s2[,1:n.snp]%*%Gamma_est
# prs_m = genotype_s2[,1:n.snp]%*%alpha_est
# model = lm(prs_y~prs_m)
# summary(model)
# 
# 
