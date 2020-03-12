#test the mediation model part


Regression = function(Y,M,A,C){
  n = length(Y)
  model1 = lm(Y~A+M+C-1)
  xf = cbind(A,M,C)
  sum.model1 = summary(model1)
  coef_temp = coef(sum.model1)
  Gamma = coef_temp[1,1]
  var_Gamma = coef_temp[1,2]^2
  model2 = lm(Y~A+C-1)
  xr = cbind(A,C)
  coef_temp2 = coef(summary(model2))
  gamma = coef_temp2[1,1]
  var_gamma = coef_temp2[1,2]^2
  p = solve(t(xf)%*%xf)%*%t(xf)%*%xr%*%solve(t(xr)%*%xr)
  cov_est = var_Gamma+var_gamma-2*p[1,1]*sum.model1$sigma^2
  return(c(Gamma,var_Gamma,
           gamma,var_gamma,
           gamma-Gamma,
           cov_est))
  
  
}

times = 2000
Gamma_est <- rep(0,times)
Gamma_var <- rep(0,times)
gamma_est <- rep(0,times)
gamma_var <- rep(0,times)
ratio_est <- rep(0,times)
ratio_var <- rep(0,times)
set.seed(666)
for(i in 1:times){
  n = 15000
  beta_a = 0.5
  beta_M = 0.3
  beta_c = 0.05
  alpha_a = 0.3
  alpha_c = 0.03
  sigma_m = 1
  sigma_y = 1
  A = rnorm(n)
  C = rnorm(n)
  M = alpha_a*A+alpha_c*C+rnorm(n,sd=sqrt(sigma_m))
  Y = beta_a*A+beta_M*M+beta_c*C+rnorm(n,sd=sqrt(sigma_y))
  temp_est = Regression(Y,M,A,C)
  Gamma_est[i] <- temp_est[1]
  Gamma_var[i] <- temp_est[2]
  gamma_est[i] <- temp_est[3]
  gamma_var[i] <- temp_est[4]
  ratio_est[i] <- temp_est[5]
  ratio_var[i] <- temp_est[6]
  
  
}









#test the mediation model part

Regression1 = function(Y,M,A){
  n = length(Y)
  model1 = lm(Y~A-1)
  sum.model1 = summary(model1)
  coef_temp = coef(sum.model1)
  Gamma = coef_temp[1,1]
  var_Gamma = coef_temp[1,2]^2
  model2 = lm(M~A-1)
  coef_temp2 = coef(summary(model2))
  gamma = coef_temp2[1,1]
  var_gamma = coef_temp2[1,2]^2
  # p = solve(t(xf)%*%xf)%*%t(xf)%*%xr%*%solve(t(xr)%*%xr)
  cov_est = solve(t(A)%*%A)%*%Gamma/gamma*var_gamma*(n-1)
  return(c(Gamma,var_Gamma,
           gamma,var_gamma,cov_est))
  
  
}


Regression2 = function(Y,M,A,A2){
  n = length(Y)
  model1 = lm(Y~A-1)
  sum.model1 = summary(model1)
  coef_temp = coef(sum.model1)
  Gamma = coef_temp[1,1]
  var_Gamma = coef_temp[1,2]^2
  model2 = lm(M~A2-1)
  coef_temp2 = coef(summary(model2))
  gamma = coef_temp2[1,1]
  var_gamma = coef_temp2[1,2]^2
 # p = solve(t(xf)%*%xf)%*%t(xf)%*%xr%*%solve(t(xr)%*%xr)
  cov_est = solve(t(A)%*%A)%*%Gamma/gamma*var_gamma*(n-1)
  return(c(Gamma,var_Gamma,
           gamma,var_gamma,cov_est))
  
  
}

times = 2000
Gamma_est <- rep(0,times)
Gamma_var <- rep(0,times)
gamma_est <- rep(0,times)
gamma_var <- rep(0,times)
ratio_est <- rep(0,times)
ratio_var <- rep(0,times)
set.seed(666)
for(i in 1:times){
  n = 1500
  beta_a = 0.5
  beta_M = 10
  #beta_c = 0.05
  alpha_a = 0.3
  #alpha_c = 0.03
  sigma_m = 1
  sigma_y = 1
  A = rnorm(n)
  C = rnorm(n)
  M = alpha_a*A+rnorm(n,sd=sqrt(sigma_m))
  Y = beta_M*M+rnorm(n,sd=sqrt(sigma_y))
  temp_est = Regression1(Y,M,A)
  # A2 = rnorm(n)
  # M2 = alpha_a*A2+rnorm(n,sd=sqrt(sigma_m))
  # temp_est = Regression2(Y,M2,A,A2)
  Gamma_est[i] <- temp_est[1]
  Gamma_var[i] <- temp_est[2]
  gamma_est[i] <- temp_est[3]
  gamma_var[i] <- temp_est[4]
  cov_est[i] <- temp_est[5]
  
  
}




Regression3 = function(Y,M,A){
  n = length(Y)
  model1 = lm(Y~A-1)
  sum.model1 = summary(model1)
  coef_temp = coef(sum.model1)
  Gamma = coef_temp[1,1]
  var_Gamma = coef_temp[1,2]^2
  model2 = lm(M~A-1)
  coef_temp2 = coef(summary(model2))
  gamma = coef_temp2[1,1]
  var_gamma = coef_temp2[1,2]^2
  n.simu = 100000
  Z1 = rnorm(n.simu,mean = Gamma,sd = sqrt(var_Gamma))
  Z2 = rnorm(n.simu,mean = gamma,sd =sqrt(var_gamma))
  distribution = Z1*Z2
  q_temp = quantile(distribution,c(0.025,0.975))
  cover = ifelse(Gamma*gamma>=q_temp[1]&
                 Gamma*gamma<=q_temp[2],1,0)
  return(c(Gamma,var_Gamma,
           gamma,var_gamma,
           cover))
  
  
}

times = 2000
Gamma_est <- rep(0,times)
Gamma_var <- rep(0,times)
gamma_est <- rep(0,times)
gamma_var <- rep(0,times)
cover <- rep(0,times)

set.seed(666)
for(i in 1:times){
  n = 15000
  beta_a = 1
  alpha_a = 2
  sigma_m = 1
  sigma_y = 1
  A = rnorm(n)
  C = rnorm(n)
  M = alpha_a*A+rnorm(n,sd=sqrt(sigma_m))
  Y = beta_a*A+rnorm(n,sd=sqrt(sigma_y))
  temp_est = Regression3(Y,M,A)
  Gamma_est[i] <- temp_est[1]
  Gamma_var[i] <- temp_est[2]
  gamma_est[i] <- temp_est[3]
  gamma_var[i] <- temp_est[4]
  cover[i] <- temp_est[5]
  
}
mean(cover)






