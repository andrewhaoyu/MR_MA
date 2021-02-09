#account for inner data correlation
IVW_c_inner = function(Gamma,var_Gamma,gamma,var_gamma){
  p <- length(Gamma)
  raio_vec = rep(0,p)
  ratio_var_vec = rep(0,p)
  
  raio_vec = Gamma/gamma
  ratio_var_vec =  var_Gamma/gamma^2+var_gamma*Gamma^2/gamma^4-2*Gamma^2*var_gamma/gamma^4
  
  
  Meta_result = Meta(raio_vec,ratio_var_vec)
  ratio_ivw =   Meta_result[1]
  ratio_ivw_var = Meta_result[2]
  coef_low = ratio_ivw-1.96*sqrt(ratio_ivw_var)
  coef_high = ratio_ivw+1.96*sqrt(ratio_ivw_var)
  cover = ifelse((beta_M>=coef_low&
                   beta_M<=coef_high),1,0)
  p_value = 2*pnorm(-abs(ratio_ivw/sqrt(ratio_ivw_var)),lower.tail = T)
  
  return(c(ratio_ivw,ratio_ivw_var,
           coef_low,coef_high,p_value,cover))
}

IVW_c = function(Gamma,var_Gamma,gamma,var_gamma){
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
  cover = ifelse((beta_M>=coef_low&
                   beta_M<=coef_high),1,0)
  p_value = 2*pnorm(-abs(ratio_ivw/sqrt(ratio_ivw_var)),lower.tail = T)
  
  return(c(ratio_ivw,ratio_ivw_var,
           coef_low,coef_high,p_value,cover))
}
Meta = function(coef_vec,var_vec){
  meta_var = (sum(1/var_vec))^-1
  meta_coef = meta_var*sum(coef_vec/var_vec)
  return(c(meta_coef,meta_var))
}
# MRLR <- function(Gamma,var_Gamma,gamma,var_gamma){
#   K <- length(Gamma)
#   keep.ind <- c(1:K)
#   
#   #first step
#   model1 = lm(Gamma~gamma-1)
#   coef_est = coefficients(model1)
#   W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
#   
#   coef_best = sum(Gamma*gamma*W_vec)/sum(gamma^2*W_vec)
#   sigma_est  = sum((Gamma-coef_est*gamma)^2)/(K-1)
#   
#   W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
#   xwx_iv = 1/sum(gamma^2*W_vec)
#   
#   var_coef_est = sigma_est*xwx_iv*t(gamma)%*%diag(W_vec)%*%diag(W_vec)%*%gamma*xwx_iv
#   
#   coef_low <- coef_est+qt(0.025,(K-1))*sqrt(var_coef_est)
#   coef_high <- coef_est+qt(0.975,(K-1))*sqrt(var_coef_est)
#   #coef_low_update <- confint(model1,level=0.95)[1]
#   #coef_high_update <- confint(model1,level=0.95)[2]
#   # cover <- ifelse((beta_M>=coef_low&
#   #                    beta_M<=coef_high),1,0)
#   
#   return(c(coef_est,coef_low,coef_high))
# }


args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
#generate phenotypes for the two-sample MR analysis
n.rep = 2
n.snp.vec = c(100,1000,5000)
beta_est_result = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_cover = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_prs = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_prs_cover = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_prs_cover_sw = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_inner = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_inner_cover = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_inner_cover_sw = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_inner_prs = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_inner_prs_cover = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_inner_prs_cover_sw = matrix(0,n.rep,length(n.snp.vec))
beta_est_IVW = matrix(0,n.rep,length(n.snp.vec))
beta_est_IVW_cover = matrix(0,n.rep,length(n.snp.vec))
beta_est_IVW_inner = matrix(0,n.rep,length(n.snp.vec))
beta_est_IVW_inner_cover = matrix(0,n.rep,length(n.snp.vec))
beta_est_IVW_inner = matrix(0,n.rep,length(n.snp.vec))
beta_est_IVW_inner_cover = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_summary = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_summary_cover = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_summary_cover_sw = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_inner_summary = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_inner_summary_cover = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_inner_summary_cover_sw = matrix(0,n.rep,length(n.snp.vec))
#beta_est_MRLR = matrix(0,n.rep,length(n.snp.vec))
#beta_est_MRLR_inner = matrix(0,n.rep,length(n.snp.vec))
set.seed(i1)
for(m in 1:length(n.snp.vec)){
  #load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata")
  load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata")
  n.sub <- nrow(genotype_s)
  n.snp = n.snp.vec[m]
  genotype_s = genotype_s[,1:n.snp]
  
  beta_M = 0.4
  h2  = 0.4
  #sigma_m = 1
  #sigma_u = alpha_U^2*var_U
  #sigma_e = 1-sigma_G-sigma_u
  sigma_m = 1-h2
  #generate M phenotypes
  
  alpha_G = rnorm(n.snp,sd=sqrt(h2/n.snp))
  #alpha_G[1] =0.1
  
  M_mat <- matrix(0,n.sub,n.rep)
  #G_value = genotype_s%*%alpha_G + alpha_U*U
  G_value = genotype_s%*%alpha_G 
  for(j in 1:n.rep){
    print(j)
    M_mat[,j] = G_value+rnorm(n.sub,sd = sqrt(sigma_m))
    #M_mat[,j] = G_value
    #+rnorm(n.sub,sd = sqrt(sigma_e))
  }
  # save(M_mat,file = "/data/zhangh24/MR_MA/result/simulation/prs/M_mat.rdata")
  
  
  #generate Y phenotypes
  #load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_Y_500k.rdata")
  load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_Y.rdata")
  genotype_s2 = genotype_s2[,1:n.snp]
  n.sub <- nrow(genotype_s2)
  n.rep = 2
  M_mat_inner = matrix(0,n.sub,n.rep)
  Y_mat <- matrix(0,n.sub,n.rep)
  G_value2 =  genotype_s2%*%alpha_G
  #+ alpha_U*U2
  #G_value2 =  genotype_s2%*%alpha_G 
  #+ alpha_U*U2
  #sigma_ey = sigma_G*beta_M^2/0.2-beta_M^2-beta_U^2*var_U
  rho = 0.3
  sigma_y = 1-beta_M^2
  sigma_ym = sqrt(sigma_y*sigma_m)*rho
  Sigma = matrix(c(sigma_y,sigma_ym,sigma_ym,sigma_m),2,2)
  #set.seed(i1)
  library(MASS)
  for(j in 1:n.rep){
    print(j)
    error = mvrnorm(n.sub,mu = c(0,0),Sigma = Sigma)
    M = G_value2+error[,2]
    Y_mat[,j] = beta_M*M+error[,1]
    M_mat_inner[,j] = M
    # sigma_e = sigma_m
    # sigma_ey = sigma_y
    # 
    # M =G_value2 + rnorm(n.sub,sd = sqrt(sigma_e))  
    # M_mat_inner[,j] = M
    # Y_mat[,j] = beta_M*M+rnorm(n.sub,sd=sqrt(sigma_ey))  
    # 
    # M =G_value2 + rnorm(n.sub,sd = sqrt(sigma_e))  
    # M_mat_inner[,j] = M
    # #M =G_value2 
    # #+ rnorm(n.sub,sd = sqrt(sigma_e))  
    # Y_mat[,j] = beta_M*M+rnorm(n.sub,sd=sqrt(sigma_y))  
    # 
    #M =G_value2 
    #+ rnorm(n.sub,sd = sqrt(sigma_e))  
    #Y_mat[,j] = beta_M*M
    #+rnorm(n.sub,sd=sqrt(sigma_ey))  
    #Y_mat[,j] = beta_M*M+beta_U*U2+rnorm(n.sub,sd=sqrt(sigma_ey))  
  }
  
  library(RcppArmadillo)
  FitLinearmodel <- function(y,x){
    model <- fastLm(X=x,y=y)
    if(is.na(coef(model)[2])){
      result <- c(0,1,1)
    }else{
      result <- coef(summary(model))[2,c(1,2,4)]  
    }
    return(result)
  }
  
  n.rep = 2
  
  beta_est = rep(0,n.snp)
  alpha_est = rep(0,n.snp)
  alpha_sd = rep(0,n.snp)
  alpha_p = rep(0,n.snp)
  alpha_est_inner = rep(0,n.snp)
  alpha_sd_inner = rep(0,n.snp)
  alpha_p_inner = rep(0,n.snp)
  Gamma_est = rep(0,n.snp)
  Gamma_sd = rep(0,n.snp)
  Gamma_p = rep(0,n.snp)
  n.train = 100000
  #n.train = 500000
  M_mat_train = M_mat[(1:n.train),]
  genotype_m_train = genotype_s[(1:n.train),]
  
  
  beta_est = rep(0,n.rep)
  
  for(l in 1:n.rep){
    for(k in 1:n.snp){
      print(k)
      #model1 <- lm(M_mat_train[,2]~genotype_m_train[,k])
      #alpha_est[k] = coef(summary(model1))[2,1]
      alpha_result = FitLinearmodel(M_mat_train[,l],cbind(1,genotype_m_train[,k]))
      alpha_est[k] = alpha_result[1]
      alpha_sd[k] = alpha_result[2]
      alpha_p[k] = alpha_result[3]
      alpha_result_inner = FitLinearmodel(M_mat_inner[,l],cbind(1,genotype_s2[,k]))
      alpha_est_inner[k] = alpha_result_inner[1]
      alpha_sd_inner[k] = alpha_result_inner[2]
      alpha_p_inner[k] = alpha_result_inner[3]
      #model2 <- lm(Y_mat[,2]~genotype_s2[,k])
      #Gamma_est[k] = coef(summary(model2))[2,1]
      Gamma_result = FitLinearmodel(Y_mat[,l],cbind(1,genotype_s2[,k]))
      Gamma_est[k] = Gamma_result[1]
      Gamma_sd[k] = Gamma_result[2]
      Gamma_p[k] = Gamma_result[3]
      #beta_est[k] = Gamma_est/alpha_est
    }
    
    genotype_m_test = genotype_s2
    
    #prs method
    prs_m_mat <- genotype_m_test%*%alpha_est
    prs_m_mat_inner = genotype_m_test%*%alpha_est_inner
    prs_y_mat <- genotype_m_test%*%Gamma_est
    #two sample MR-PRS using individual level data
   
    #model1 = lm(prs_y_mat~prs_m_mat)
    #coefficients(model1)
    #crossprod(prs_y_mat,prs_m_mat)/crossprod(prs_m_mat)
    #model2 = lm(Y_mat[,l]~prs_m_mat)
    #coefficients(model2)
    model = lm(Y_mat[,l]~prs_m_mat)
    prs_m_train = genotype_m_train%*%alpha_est
    model_m = lm(M_mat_train[,l]~prs_m_train)
    sigma_m_est = summary(model_m)$sigma^2
    #sigma_m_est = sigma_m
    F = crossprod(G_value2)/n.snp/sigma_m_est
    
    beta_est_result[l,m]= coefficients(summary(model))[2,1]*(F+1)/F
    beta_est = as.numeric(beta_est_result[l,m])
    sigma_y_est = var(Y_mat[,l]-M_mat_train[,l]*beta_est)
    beta_temp = beta_est
    
    beta_est_var = sigma_y_est/crossprod(prs_m_mat)
    beta_est_result_low <- beta_est-1.96*sqrt(beta_est_var)
    beta_est_result_high <- beta_est+1.96*sqrt(beta_est_var)
    beta_est_result_cover[l,m] <- ifelse((beta_M>=beta_est_result_low)&
                                      (beta_M<=beta_est_result_high),1,0)
    
    
    model = lm(prs_y_mat~prs_m_mat)
    beta_est_result_prs[l,m]= coefficients(summary(model))[2,1]*(F+1)/F
    beta_est = as.numeric(beta_est_result_prs[l,m])
    beta_est_var = var(Y_mat[,l]-M_mat_train[,l]*beta_est)/crossprod(prs_m_mat)
    sigma_my_est = as.numeric(cov(Y_mat[,l]-beta_est*M_mat_inner[,l],M_mat_inner[,l]-prs_m_mat_inner))
    beta_est_result_low <- beta_est-1.96*sqrt(beta_est_var)
    beta_est_result_high <- beta_est+1.96*sqrt(beta_est_var)
    beta_est_result_prs_cover[l,m] <- ifelse((beta_M>=beta_est_result_low)&
                                           (beta_M<=beta_est_result_high),1,0)
    beta_est_var = (sigma_y_est+beta_est*sigma_my_est)/crossprod(prs_m_mat)
    beta_est_result_low <- beta_est-1.96*sqrt(beta_est_var)
    beta_est_result_high <- beta_est+1.96*sqrt(beta_est_var)
    beta_est_result_prs_cover[l,m] <- ifelse((beta_M>=beta_est_result_low)&
                                              (beta_M<=beta_est_result_high),1,0)
    
    #two-sample MR using summary level data
    sigma_m_est = 1-sum(alpha_est^2-alpha_sd^2)
    #sigma_m_est = 1-sum(alpha_est^2)
    #(crossprod(alpha_est,Gamma_est)/crossprod(alpha_est))
    # prs_m = genotype_m_test%*%alpha_est
    # prs_y = genotype_m_test%*%Gamma_est
    # crossprod(prs_y,prs_m)/crossprod(prs_m,prs_m)
    #sigma_m_est = sigma_m
    N <- nrow(prs_y_mat)
    F = N*sum(alpha_est^2-alpha_sd^2)/sigma_m_est/n.snp
    #F = N*sum(alpha_est^2)/sigma_m_est/n.snp
    beta_est_result_summary[l,m] <- (crossprod(alpha_est,Gamma_est)/crossprod(alpha_est))*(F+1)/F
    beta_est = as.numeric(beta_est_result_summary[l,m])
    
    beta_est_var = (1-beta_est^2)/(N*sum(alpha_est^2-alpha_sd^2)) 
    beta_est_result_summary_low <- beta_est-1.96*sqrt(beta_est_var)
    beta_est_result_summary_high <- beta_est+1.96*sqrt(beta_est_var)
    beta_est_result_summary_cover[l,m] <- ifelse((beta_M>=beta_est_result_summary_low)&
                                           (beta_M<=beta_est_result_summary_high),1,0)
    
    
    #one sample analysis using individual level data
    prs_m_mat_inner = genotype_m_test%*%alpha_est_inner
    #model  = lm(prs_y_mat~prs_m_mat_inner)
    model  = lm(Y_mat[,l]~prs_m_mat_inner)
    model_m = lm(M_mat_inner[,l]~prs_m_mat_inner)
    sigma_m_est = summary(model_m)$sigma^2
    #sigma_m_est = sigma_m
    F = crossprod(G_value2)/n.snp/sigma_m_est
    beta_temp = as.numeric(coefficients(summary(model))[2,1])
    #sigma_my = cov(Y_mat[,l]-M_mat_inner[,l]*beta_temp,
     #              model_m$residuals)
    sigma_my = as.numeric(cov(Y_mat[,l]-beta_temp*M_mat_inner[,l],M_mat_inner[,l]-prs_m_mat_inner))
    beta_est_result_inner[l,m] = coefficients(summary(model))[2,1]-as.numeric(sigma_my/(sigma_m_est*(F+1)))
    beta_est = as.numeric(beta_est_result_inner[l,m])
    
    
    beta_est_var = var(Y_mat[,l]-M_mat_inner[,l]*beta_est)/crossprod(prs_m_mat_inner)
    beta_est_result_low <- beta_est-1.96*sqrt(beta_est_var)
    beta_est_result_high <- beta_est+1.96*sqrt(beta_est_var)
    beta_est_result_inner_cover[l,m] <- ifelse((beta_M>=beta_est_result_low)&
                                           (beta_M<=beta_est_result_high),1,0)
    model  = lm(prs_y_mat~prs_m_mat_inner)
    #model  = lm(Y_mat[,l]~prs_m_mat_inner)
    model_m = lm(M_mat_inner[,l]~prs_m_mat_inner)
    sigma_m_est = summary(model_m)$sigma^2
    #sigma_m_est = sigma_m
    F = crossprod(G_value2)/n.snp/sigma_m_est
    beta_temp = as.numeric(coefficients(summary(model))[2,1])
    #sigma_my = cov(Y_mat[,l]-M_mat_inner[,l]*beta_temp,
    #              model_m$residuals)
    sigma_my_est = as.numeric(cov(Y_mat[,l]-beta_temp*M_mat_inner[,l],M_mat_inner[,l]-prs_m_mat_inner))
    beta_est_result_inner_prs[l,m] = coefficients(summary(model))[2,1]-as.numeric(sigma_my_est/(sigma_m_est*(F+1)))
    beta_est = as.numeric(beta_est_result_inner_prs[l,m])
    sigma_y_est = var(Y_mat[,l]-M_mat_inner[,l]*beta_est)
    beta_est_var = (sigma_y_est)/crossprod(prs_m_mat_inner)
    beta_est_result_low <- beta_est-1.96*sqrt(beta_est_var)
    beta_est_result_high <- beta_est+1.96*sqrt(beta_est_var)
    beta_est_result_inner_prs_cover[l,m] <- ifelse((beta_M>=beta_est_result_low)&
                                                 (beta_M<=beta_est_result_high),1,0)
    
    beta_est_var = (sigma_y_est+sigma_my_est*beta_est)/crossprod(prs_m_mat_inner)
    beta_est_result_low <- beta_est-1.96*sqrt(beta_est_var)
    beta_est_result_high <- beta_est+1.96*sqrt(beta_est_var)
    beta_est_result_inner_prs_cover_sw[l,m] <- ifelse((beta_M>=beta_est_result_low)&
                                                     (beta_M<=beta_est_result_high),1,0)
    
    
    
    #one sample analysis using summary level data
    sigma_m_est = 1-sum(alpha_est^2)
    #sigma_m_est = sigma_m
    N <- nrow(prs_y_mat)
    #model  = lm(prs_y_mat~prs_m_mat_inner)
    beta_est_result_inner_summary[l,m] = (crossprod(alpha_est_inner,Gamma_est)/crossprod(alpha_est_inner))-sigma_my/(sigma_m_est*(F+1))
    beta_est = as.numeric(beta_est_result_inner_summary[l,m])
    #beta_est_var = (1-beta_est^2)/(N*sum(alpha_est^2-alpha_sd^2)) 
    beta_est_var = (1-beta_est^2)/(N*sum(alpha_est^2)) 
    beta_est_result_low <- beta_est-1.96*sqrt(beta_est_var)
    beta_est_result_high <- beta_est+1.96*sqrt(beta_est_var)
    beta_est_result_inner_summary_cover[l,m] <- ifelse((beta_M>=beta_est_result_low)&
                                                 (beta_M<=beta_est_result_high),1,0)
    
    
    #IVW method
    idx <- which(alpha_p<=5E-08)
    
    beta_est_IVW[l,m] = IVW_c(Gamma_est[idx],Gamma_sd[idx]^2,alpha_est[idx],alpha_sd[idx]^2)[1]
    beta_est_IVW_cover[l,m] = IVW_c(Gamma_est[idx],Gamma_sd[idx]^2,alpha_est[idx],alpha_sd[idx]^2)[6]
    #beta_est_MRLR[l,m] = MRLR(Gamma_est[idx],Gamma_sd[idx]^2,alpha_est[idx],alpha_sd[idx]^2)[1]
    
    idx <- which(alpha_p_inner<=5E-08)
    beta_est_IVW_inner[l,m] = IVW_c_inner(Gamma_est[idx],Gamma_sd[idx]^2,alpha_est_inner[idx],alpha_sd_inner[idx]^2)[1]
    beta_est_IVW_inner_cover[l,m] = IVW_c_inner(Gamma_est[idx],Gamma_sd[idx]^2,alpha_est_inner[idx],alpha_sd_inner[idx]^2)[6]
    #beta_est_MRLR_inner[l,m] = MRLR(Gamma_est[idx],Gamma_sd[idx]^2,alpha_est_inner[idx],alpha_sd_inner[idx]^2)[1]
    
  }
  
}
result= list(beta_est_result,
             beta_est_result_cover,
             beta_est_result_summary,
             beta_est_result_summary_cover,
             beta_est_result_inner,
             beta_est_result_inner_cover,
             beta_est_result_inner_summary,
             beta_est_result_inner_summary_cover,
             beta_est_IVW,
             beta_est_IVW_cover,
             beta_est_IVW_inner,
             beta_est_IVW_inner_cover,
             beta_est_result_prs,
             beta_est_result_prs_cover,
             beta_est_result_inner_prs,
             beta_est_result_inner_prs_cover,
             beta_est_result_prs_cover_sw,
             beta_est_result_inner_prs_cover_sw)
save(result,file  = paste0("/data/zhangh24/MR_MA/result/simulation/prs/beta_test_result_rho_",i1))
# confint(model)
# 
# library(ggplot2)
# data = data.frame(alpha= alpha_est[idx],Gamma = Gamma_est[idx])
# ggplot(data,aes(alpha,Gamma))+ 
#   geom_point()+
#   geom_smooth(method="lm")+
#   geom_abline(slope=0.15,
#               intercept=0,
#               col = "red")+
#   theme_Publication()+
#   xlab("alpha estiamte")+
#   ylab("Gamma estimate")+
#   ggtitle("Two sample MR")
# 
# idx <- which(alpha_p_inner<=5E-08)
# data = data.frame(alpha= alpha_est_inner[idx],Gamma = Gamma_est[idx])
# 
# ggplot(data,aes(alpha,Gamma))+ 
#   geom_point()+
#   geom_smooth(method="lm")+
#   geom_abline(slope=0.15,
#               intercept=0,
#               col = "red")+
#   theme_Publication()+
#   xlab("alpha estiamte")+
#   ylab("Gamma estimate")+
#   ggtitle("One sample MR")
# 
# 
# idx <- which(alpha_p<=5E-08)
# plot(alpha_est[idx],Gamma_est[idx])
# abline(a=0,b=0.15)
# 
# 
# prs_y = genotype_s2[,1:n.snp]%*%Gamma_est
# prs_m = genotype_s2[,1:n.snp]%*%alpha_est
# model = lm(prs_y~prs_m)
# summary(model)
# save(Y_mat,file = "/data/zhangh24/MR_MA/result/simulation/prs/Y_mat.rdata")

