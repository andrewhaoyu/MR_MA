MRLR <- function(Gamma,se_Gamma,gamma,se_gamma,R,L){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  #first step
  model1 = lm(Gamma~gamma-1)
  coef_est = coefficients(model1)
  
  #
  
  
  #W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  W_inv = GetWMat(var_Gamma,var_gamma,
                  se_Gamma,se_gamma,coef_est,R)
  coef_best = (t(gamma)%*%W_inv%*%Gamma)/(t(gamma)%*%W_inv%*%gamma)
  
  coef_vec <- seq(from = coef_best-0.1,to = coef_best+0.1,by  = 0.001)
  
  quad_vec = rep(0,length(coef_vec))
  for(k in 1:length(coef_vec)){
    coef_temp = coef_vec[k]
    W_inv = GetWMat(var_Gamma,var_gamma,
                    se_Gamma,se_gamma,coef_temp,R)
    
    quad_vec[k] = t(Gamma-coef_temp*gamma)%*%W_inv%*%(Gamma-coef_temp*gamma)
  }
  
  coef_best = coef_vec[which.min(quad_vec)]
  return(coef_best)
}
