N <- 1000
library(MASS)
Sigma = matrix(c(1,0.8,0.8,1),2)
X <- mvrnorm(N,mu = c(0,0),Sigma = Sigma)
beta = 0.15
sigma_y = 1- beta^2
sigma_m = 0.4
rho = 0.3
sigma_my = rho*(sigma_y)*sigma_m
Sigma_e = matrix(c(sigma_y,sigma_my,sigma_my,sigma_m),2,2)
error = mvrnorm(N,mu = c(0,0),Sigma = Sigma)
Y = lm()
