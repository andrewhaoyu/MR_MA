GetWMat <- function(var_Gamma,var_gamma,
                    se_Gamma,se_gamma,coef_est,R){
  K = length(Gamma)
  W  = matrix(0,K,K)
  diag(W) = var_Gamma+coef_est^2*var_gamma
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      W[j,i] = W[i,j] = se_Gamma[i]*se_Gamma[j]*R[i,j]+coef_est^2*se_gamma[i]*se_gamma[j]
    }
  }
  W_inv = solve(W)
  return(W_inv)
}


MRLR <- function(Gamma,var_Gamma,gamma,var_gamma,R){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  se_Gamma  = sqrt(var_Gamma)
  se_gamma = sqrt(var_gamma)
  #first step
  model1 = lm(Gamma~gamma-1)
  coef_est = coefficients(model1)
  
 
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

TestPleotropic <- function(Gamma,var_Gamma,gamma,var_gamma,coef_best){
  beta_plug = coef_best
  
  Gamma_update = Gamma
  var_Gamma_update = var_Gamma
  gamma_update = gamma
  var_gamma_update = var_gamma
  p_test_value = pchisq((Gamma_update-beta_plug*gamma_update)^2/(var_Gamma_update+beta_plug^2*var_gamma_update),1,lower.tail = F)
  p_adjust_test_value = p.adjust(p_test_value,method="none")  
  return(p_adjust_test_value)
}


MRWeight <- function(Gamma,var_Gamma,gamma,var_gamma,R){
  coef_best = MRLR(Gamma,var_Gamma,gamma,var_gamma,R)
  #test pleotropic
  K <- length(Gamma)
  all.id <- c(1:K)
  #id of SNPs that are removed
  out.id <- NULL
  #id of SNPs that are kept
  keep.id <- all.id
  
  p_adjust_test_value <- TestPleotropic(Gamma,var_Gamma,gamma,var_gamma,coef_best)
  cut.off = 0.01
  if(min(p_adjust_test_value)<cut.off){
    while(min(p_adjust_test_value)<cut.off){
      
      out.id <- c(out.id,keep.id[which.min(p_adjust_test_value)])
      keep.id <- setdiff(all.id,out.id)
      Gamma.keep = Gamma[keep.id]
      var_Gamma.keep = var_Gamma[keep.id]
      gamma.keep = gamma[keep.id]
      var_gamma.keep = var_gamma[keep.id]
      coef_best = MRLR(Gamma.keep,
                       var_Gamma.keep,
                       gamma.keep,
                       var_gamma.keep,
                       R)
      p_adjust_test_value <- TestPleotropic(Gamma.keep,
                                            var_Gamma.keep,
                                            gamma.keep,
                                            var_gamma.keep,
                                            coef_best)
    }
    
    
  }else{
    Gamma.keep = Gamma[keep.id]
    var_Gamma.keep = var_Gamma[keep.id]
    gamma.keep = gamma[keep.id]
    var_gamma.keep = var_gamma[keep.id]
  }
  coef_est = coef_best
  result.temp <- MRLRVariance(Gamma.keep,
                              var_Gamma.keep,
                              gamma.keep,
                              var_gamma.keep,
                              coef_best,R)
  coef_var = result.temp[1]
  coef_low = result.temp[2]
  coef_high = result.temp[3]
  #coef_low_update <- confint(model1,level=0.95)[1]
  #coef_high_update <- confint(model1,level=0.95)[2]
  # cover <- ifelse((beta_M>=coef_low&
  #                    beta_M<=coef_high),1,0)
  # 
  return(list(coef_est,coef_low,coef_high,sqrt(coef_var),keep.id, out.id))
}

MRLRVariance <- function(Gamma,var_Gamma,gamma,var_gamma,coef_est,R){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  
  sigma_est  = sum((Gamma-coef_est*gamma)^2)/(K-1)
  W_inv = GetWMat(var_Gamma,var_gamma,
                  se_Gamma,se_gamma,coef_est,R)
  
  
  xwx_iv = solve(t(gamma)%*%W_inv%*%gamma)
  var_coef_est = sigma_est*xwx_iv*t(gamma)%*%W_inv%*%W_inv%*%gamma*xwx_iv
  coef_low <- coef_est+qt(0.025,(K-1))*sqrt(var_coef_est)
  coef_high <- coef_est+qt(0.975,(K-1))*sqrt(var_coef_est)
  return(c(var_coef_est,
           coef_low,
           coef_high))
}






#sqrt(MRLRVariance(Gamma,var_Gamma,gamma,var_gamma,coef_best,R)[1])
