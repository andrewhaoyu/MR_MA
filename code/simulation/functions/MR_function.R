Meta = function(coef_vec,var_vec){
  meta_var = (sum(1/var_vec))^-1
  meta_coef = meta_var*sum(coef_vec/var_vec)
  return(c(meta_coef,meta_var))
}
MetaWeight = function(coef_vec,var_vec){
  meta_var = (sum(1/var_vec))^-1
  weight = meta_var/var_vec
  return(weight)
}

MRLR <- function(Gamma,var_Gamma,gamma,var_gamma){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  
  #first step
  model1 = lm(Gamma~gamma-1)
  coef_est = coefficients(model1)
  W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  
  coef_best = sum(Gamma*gamma*W_vec)/sum(gamma^2*W_vec)
  
  coef_vec <- seq(coef_best-0.1,coef_best+0.1,0.001)
  
  quad_vec = rep(0,length(coef_vec))
  for(i in 1:length(coef_vec)){
    coef_temp = coef_vec[i]
    W_vec = 1/(var_Gamma+coef_temp^2*var_gamma)
    quad_vec[i] = sum((Gamma-coef_temp*gamma)^2*W_vec)
  }
  
  coef_best = coef_vec[which.min(quad_vec)]
  return(coef_best)
}
MRLRVariance <- function(Gamma,var_Gamma,gamma,var_gamma,coef_est){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  
  sigma_est  = sum((Gamma-coef_est*gamma)^2)/(K-1)
  
  W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  xwx_iv = 1/sum(gamma^2*W_vec)
  var_coef_est = sigma_est*xwx_iv*t(gamma)%*%diag(W_vec)%*%diag(W_vec)%*%gamma*xwx_iv
  coef_low <- coef_est+qt(0.025,(K-1))*sqrt(var_coef_est)
  coef_high <- coef_est+qt(0.975,(K-1))*sqrt(var_coef_est)
  return(c(var_coef_est,
           coef_low,
           coef_high))
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







MRWeight <- function(Gamma,var_Gamma,gamma,var_gamma){
  coef_best = MRLR(Gamma,var_Gamma,gamma,var_gamma)
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
                       var_gamma.keep)
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
                              coef_best)
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


IVW_c <- function(Gamma,var_Gamma,gamma,var_gamma){
  p <- length(Gamma)
  raio_vec = rep(0,p)
  ratio_var_vec = rep(0,p)
  
  raio_vec = Gamma/gamma
  ratio_var_vec =  var_Gamma/gamma^2+var_gamma*Gamma^2/gamma^4
  
  
  Meta_result = Meta(raio_vec,ratio_var_vec)
  ratio_ivw =   Meta_result[1]
  ratio_ivw_var = Meta_result[2]
  coef_low = ratio_ivw-1.96*sqrt(ratio_ivw_var)
  coef_high = ratio_ivw+1.96*sqrt(ratio_ivw_var)
  
  
  return(c(ratio_ivw,
           coef_low,coef_high,sqrt(ratio_ivw_var)))
}

