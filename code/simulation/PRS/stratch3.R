#simulation setting
#three SNPs
#two SNPs are independent causal SNPs
#one is correlated with both of the two SNPs
#first estimate are using two true causal SNPs
#send estimate is using the tagging SNP
n.rep = 1000
beta_M_est1 = rep(0,n.rep)
beta_M_est2 = rep(0,n.rep)
alpha_est = rep(0,n.rep)
library(MASS)

N = 10000
for(i_rep in 1:n.rep){
  print(i_rep)
  Sigma_G = matrix(c(1,0.5,0,
                     0.5,1,0.8,
                     0,0.8,1),3,3)
  G = mvrnorm(N,mu = rep(0,3),Sigma =Sigma_G)
  beta_M = 0.3
  G1 = G[,1]
  G2 = G[,2]
  G3 = G[,3]
  alpha = c(0.1,0.1)
  sigma_e  = 1- sum(alpha^2)
  sigma_y = 1- beta_M^2
  rho = 0.3
  sigma_ey = sqrt(sigma_e*sigma_y)*rho
  Sigma_e = matrix(c(sigma_e,sigma_ey,sigma_ey,sigma_y),2,2)
  error = mvrnorm(N,mu = rep(0,2),Sigma = Sigma_e)
  M = cbind(G1,G3)%*%alpha+ error[,2]
  Y = beta_M*M+error[,1]
  model1 = lm(M~cbind(G1,G3))
  model2 = lm(Y~predict(model1))
  summary(model2)
  beta_M_est1[i_rep] = coefficients(model2)[2]
  beta_M_est2[i_rep] = crossprod(G2,M)^-1*crossprod(G2,Y)
  alpha_est[i_rep] <- coefficients(lm(M~G2))[2]
}
Getestimate <- function(est){
  return(paste0(sprintf("%.3f",round(mean(est),3))," [",
                sprintf("%.3f",round(quantile(est,0.025),3)),
                ", ",
                sprintf("%.3f",round(quantile(est,0.975),3)),
                "]"))
  
  
}

Getestimate(beta_M_est1)
Getestimate(beta_M_est2)





#two SNPs
#onw SNPs is causal SNP
#the other one is correlated with both of the causal SNP


n.rep = 1000
beta_M_est1 = rep(0,n.rep)
beta_M_est2 = rep(0,n.rep)
alpha_est = rep(0,n.rep)
library(MASS)

N = 10000
for(i_rep in 1:n.rep){
  print(i_rep)
  Sigma_G = matrix(c(1,0.5,
                     0.5,1),2,2)
  G = mvrnorm(N,mu = rep(0,2),Sigma =Sigma_G)
  beta_M = 0.3
  G1 = G[,1]
  G2 = G[,2]
  alpha = c(0.1)
  sigma_e  = 1- sum(alpha^2)
  sigma_y = 1- beta_M^2
  rho = 0.3
  sigma_ey = sqrt(sigma_e*sigma_y)*rho
  Sigma_e = matrix(c(sigma_e,sigma_ey,sigma_ey,sigma_y),2,2)
  error = mvrnorm(N,mu = rep(0,2),Sigma = Sigma_e)
  M = cbind(G1)%*%alpha+ error[,2]
  Y = beta_M*M+error[,1]
  model1 = lm(M~cbind(G1))
  model2 = lm(Y~predict(model1))
  summary(model2)
  beta_M_est1[i_rep] = coefficients(model2)[2]
  beta_M_est2[i_rep] = crossprod(G2,M)^-1*crossprod(G2,Y)
  alpha_est[i_rep] <- coefficients(lm(M~G2))[2]
}

Getestimate(beta_M_est1)
Getestimate(beta_M_est2)
