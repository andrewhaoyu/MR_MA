
GetWMat <- function(var_Gamma,var_alpha,
                    se_Gamma,se_alpha,coef_est,R){
  K = length(Gamma)
  var_U  = matrix(0,K,K)
  diag(var_U) = var_Gamma+as.numeric(coef_est)^2*var_alpha
  for(i in 1:(K-1)){
    for(j in (i+1):K){
      var_U[j,i]  = var_U[i,j] = se_Gamma[i]*se_Gamma[j]*R[i,j]+coef_est^2*se_alpha[i]*se_alpha[j]*R[i,j]
    }
  }
  W = spdinv(var_U)
  return(W)
}


MRLR <- function(Gamma,var_Gamma,alpha,var_alpha,R){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  se_Gamma  = sqrt(var_Gamma)
  se_alpha = sqrt(var_alpha)
  #first step
  # model1 = lm(Gamma~gamma-1)
  # coef_est = coefficients(model1)
  # 
  
  #W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  
  loss_function <- function(coef_est){
    W = GetWMat(var_Gamma,var_alpha,
                se_Gamma,se_alpha,coef_est,R)
    loss = t(Gamma-coef_est*alpha)%*%W%*%(Gamma-coef_est*alpha)
    return(loss)
  }
  
  bound <- quantile(abs(Gamma/alpha), 0.95, na.rm = TRUE) * 
    2
  coef_best = optimize(f = loss_function,
           interval = bound * c(-1, 1), maximum = F, 
           tol = .Machine$double.eps^0.5)$minimum
  
  # coef_best = as.numeric((t(alpha)%*%W%*%Gamma)/(t(alpha)%*%W%*%alpha))
  # coef_best = beta1
  # quad1 = t(Gamma-coef_best*alpha)%*%(Gamma-coef_best*alpha) 
  # # W = GetWMat(var_Gamma,var_alpha,
  # #             se_Gamma,se_alpha,coef_best,R)
  # coef_best = beta2
  # 
  # quad2 = t(Gamma-coef_best*alpha)%*%(Gamma-coef_best*alpha) 
  #quad2 = t(Gamma-coef_est*alpha)%*%W%*%(Gamma-coef_est*alpha) 
  
  #coef_vec <- seq(from = coef_best-0.1,to = coef_best+0.1,by  = 0.001)
  
  # quad_vec = rep(0,length(coef_vec))
  # for(k in 1:length(coef_vec)){
  #   coef_temp = coef_vec[k]
  #   W = GetWMat(var_Gamma,var_alpha,
  #                   se_Gamma,se_alpha,coef_temp,R)
  #   
  #   quad_vec[k] = t(Gamma-coef_temp*alpha)%*%W%*%(Gamma-coef_temp*alpha)
  # }
  # 
  # coef_best = coef_vec[which.min(quad_vec)]
  return(coef_best)
}

TestPleotropic <- function(Gamma,var_Gamma,alpha,var_alpha,coef_best){
  beta_plug = coef_best
  
  Gamma_update = Gamma
  var_Gamma_update = var_Gamma
  alpha_update = alpha
  var_alpha_update = var_alpha
  p_test_value = pchisq((alpha_update-beta_plug*alpha_update)^2/(var_Gamma_update+beta_plug^2*var_alpha_update),1,lower.tail = F)
  p_adjust_test_value = p.adjust(p_test_value,method="none")  
  return(p_adjust_test_value)
}


MRWeight <- function(Gamma,var_Gamma,alpha,var_alpha,R){
  coef_best = MRLR(Gamma,var_Gamma,alpha,var_alpha,R)
  #test pleotropic
  K <- length(Gamma)
  all.id <- c(1:K)
  #id of SNPs that are removed
  out.id <- NULL
  #id of SNPs that are kept
  keep.id <- all.id
  
  # p_adjust_test_value <- TestPleotropic(Gamma,var_Gamma,alpha,var_alpha,coef_best)
  # cut.off = 0.01
  # if(min(p_adjust_test_value)<cut.off){
  #   while(min(p_adjust_test_value)<cut.off){
  #     
  #     out.id <- c(out.id,keep.id[which.min(p_adjust_test_value)])
  #     keep.id <- setdiff(all.id,out.id)
  #     Gamma.keep = Gamma[keep.id]
  #     var_Gamma.keep = var_Gamma[keep.id]
  #     alpha.keep = alpha[keep.id]
  #     var_alpha.keep = var_alpha[keep.id]
  #     coef_best = MRLR(Gamma.keep,
  #                      var_Gamma.keep,
  #                      alpha.keep,
  #                      var_alpha.keep,
  #                      R)
  #     p_adjust_test_value <- TestPleotropic(Gamma.keep,
  #                                           var_Gamma.keep,
  #                                           alpha.keep,
  #                                           var_alpha.keep,
  #                                           coef_best)
  #   }
  #   
  #   
  # }else{
  #   Gamma.keep = Gamma[keep.id]
  #   var_Gamma.keep = var_Gamma[keep.id]
  #   alpha.keep = alpha[keep.id]
  #   var_alpha.keep = var_alpha[keep.id]
  # }
  beta.hat = coef_best
  
  # b_exp = alpha
  # se_exp = sqrt(var_alpha)
  # se_out = sqrt(var_Gamma)
  # b_out = Gamma

  # score.var <- sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 -
  #                                                        se_out^2) * se_exp^2 + se_exp^2 * se_out^2)/(se_out^2 +
  #                                                                                                       beta.hat^2 * se_exp^2)^2)
  # I <- sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 - se_out^2) *
  #             se_exp^2)/(se_out^2 + beta.hat^2 * se_exp^2)^2)
  # 
  # coef_var = score.var/I^2
  # coef_low = coef_best - 1.96*sqrt(coef_var)
  # coef_high = coef_best + 1.96* sqrt(coef_var)
  # coef_var_backup = coef_var
  result.temp <- MRLRVariance(Gamma,
                              var_Gamma,
                              alpha,
                              var_alpha,
                              coef_best,R)
  coef_var = result.temp[1]
  coef_low = result.temp[2]
  coef_high = result.temp[3]
  # result.temp <- MRLRVariance(Gamma.keep,
  #                             var_Gamma.keep,
  #                             alpha.keep,
  #                             var_alpha.keep,
  #                             coef_best,R)
 
  #coef_low_update <- confint(model1,level=0.95)[1]
  #coef_high_update <- confint(model1,level=0.95)[2]
  # cover <- ifelse((beta_M>=coef_low&
  #                    beta_M<=coef_high),1,0)
  # 
  #return(list(coef_est,coef_low,coef_high,sqrt(coef_var),keep.id, out.id))
  return(list(coef_best,coef_low,coef_high,sqrt(coef_var)))
}

# MRLRVariance <- function(Gamma,var_Gamma,alpha,var_alpha,coef_est,R){
#   K <- length(Gamma)
#   se_Gamma = sqrt(var_Gamma)
#   se_alpha = sqrt(var_alpha)
#   
#   #sigma_est  = sum((Gamma-coef_est*alpha)^2)/(K-1)
#   W = GetWMat(var_Gamma,var_alpha,
#                   se_Gamma,se_alpha,coef_est,R)
#   
#   
#   var_coef_est = solve(t(alpha)%*%W%*%alpha)
#   #var_coef_est = sigma_est*xwx_iv*t(alpha)%*%W_inv%*%W_inv%*%alpha*xwx_iv
#   coef_low <- coef_est+qnorm(0.025,sd = sqrt(var_coef_est))
#   coef_high <- coef_est+qnorm(0.975,sd = sqrt(var_coef_est))
#   #coef_low <- coef_est+qt(0.025,(K-1))*sqrt(var_coef_est)
#   #coef_high <- coef_est+qt(0.975,(K-1))*sqrt(var_coef_est)
#   return(c(var_coef_est,
#            coef_low,
#            coef_high))
# }

MRLRVariance <- function(Gamma,var_Gamma,alpha,var_alpha,coef_est,R){
  K <- length(Gamma)
  se_Gamma = sqrt(var_Gamma)
  se_alpha = sqrt(var_alpha)
  
  sigma_est  = sum((Gamma-coef_est*alpha)^2)/(K-1)
  
  #A(beta)
  #B(beta)
  #variance = t(A(beta))^-1B(beta)A(beta)
  
  # A = sum(
  #   (alpha^2-var_alpha)/(var_Gamma+coef_est^2*var_alpha)+
  #   (2*var_alpha*coef_est*(alpha*Gamma-coef_est*(alpha^2-var_alpha)))/
  #   (var_Gamma+coef_est^2*var_alpha)^2
  #   )
  # B = sum((alpha^2-var_alpha)/(var_Gamma+coef_est^2*var_alpha))
  # 
  # 
  # A_temp = sum(
  #   (alpha^2-var_alpha)/(var_Gamma+coef_est^2*var_alpha)+
  #     (2*var_alpha*coef_est*(alpha*Gamma-sqrt(var_alpha*var_Gamma)-coef_est*(alpha^2-var_alpha)))/
  #     (var_Gamma+coef_est^2*var_alpha)^2
  # )
  # var_coef_est = B/A^2
  # 
  
  
  
  var_coef_est = solve(t(alpha)%*%W%*%alpha)
  W = GetWMat(var_Gamma,var_alpha,
              se_Gamma,se_alpha,coef_est,R)

  var_coef_est = sum((alpha^2-var_alpha)*diag(W))^-1
  #var_coef_est = sigma_est*xwx_iv*t(alpha)%*%W_inv%*%W_inv%*%alpha*xwx_iv
  coef_low <- coef_est+qnorm(0.025,sd = sqrt(var_coef_est))
  coef_high <- coef_est+qnorm(0.975,sd = sqrt(var_coef_est))
  #coef_low <- coef_est+qt(0.025,(K-1))*sqrt(var_coef_est)
  #coef_high <- coef_est+qt(0.975,(K-1))*sqrt(var_coef_est)
  return(c(var_coef_est,
           coef_low,
           coef_high))
}







IVW_c <- function(Gamma,var_Gamma,alpha,var_alpha){
  p <- length(Gamma)
  raio_vec = rep(0,p)
  ratio_var_vec = rep(0,p)
  
  raio_vec = Gamma/alpha
  ratio_var_vec =  var_Gamma/alpha^2+var_alpha*Gamma^2/alpha^4
  
  
  Meta_result = Meta(raio_vec,ratio_var_vec)
  ratio_ivw =   Meta_result[1]
  ratio_ivw_var = Meta_result[2]
  coef_low = ratio_ivw-1.96*sqrt(ratio_ivw_var)
  coef_high = ratio_ivw+1.96*sqrt(ratio_ivw_var)
  
  
  return(c(ratio_ivw,
           coef_low,coef_high,sqrt(ratio_ivw_var)))
}

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






#sqrt(MRLRVariance(Gamma,var_Gamma,alpha,var_alpha,coef_best,R)[1])

