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
  #cover = ifelse((beta_M>=coef_low&
  #                  beta_M<=coef_high),1,0)
  p_value = 2*pnorm(-abs(ratio_ivw/sqrt(ratio_ivw_var)),lower.tail = T)
  
  return(c(ratio_ivw,ratio_ivw_var,
           coef_low,coef_high,p_value))
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
  #cover = ifelse((beta_M>=coef_low&
  #                  beta_M<=coef_high),1,0)
  p_value = 2*pnorm(-abs(ratio_ivw/sqrt(ratio_ivw_var)),lower.tail = T)
  
  return(c(ratio_ivw,ratio_ivw_var,
           coef_low,coef_high,p_value))
}
Meta = function(coef_vec,var_vec){
  meta_var = (sum(1/var_vec))^-1
  meta_coef = meta_var*sum(coef_vec/var_vec)
  return(c(meta_coef,meta_var))
}
MRLR <- function(Gamma,var_Gamma,gamma,var_gamma){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  
  #first step
  model1 = lm(Gamma~gamma-1)
  coef_est = coefficients(model1)
  W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  
  coef_best = sum(Gamma*gamma*W_vec)/sum(gamma^2*W_vec)
  sigma_est  = sum((Gamma-coef_est*gamma)^2)/(K-1)
  
  W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  xwx_iv = 1/sum(gamma^2*W_vec)
  
  var_coef_est = sigma_est*xwx_iv*t(gamma)%*%diag(W_vec)%*%diag(W_vec)%*%gamma*xwx_iv
  
  coef_low <- coef_est+qt(0.025,(K-1))*sqrt(var_coef_est)
  coef_high <- coef_est+qt(0.975,(K-1))*sqrt(var_coef_est)
  #coef_low_update <- confint(model1,level=0.95)[1]
  #coef_high_update <- confint(model1,level=0.95)[2]
  # cover <- ifelse((beta_M>=coef_low&
  #                    beta_M<=coef_high),1,0)
  
  return(c(coef_est,coef_low,coef_high))
}

args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
#generate phenotypes for the two-sample MR analysis
n.rep = 2
n.snp.vec = c(100,1000,5000)
beta_est_result = matrix(0,n.rep,length(n.snp.vec))
beta_est_result_inner = matrix(0,n.rep,length(n.snp.vec))
beta_est_IVW = matrix(0,n.rep,length(n.snp.vec))
beta_est_IVW_inner = matrix(0,n.rep,length(n.snp.vec))
beta_est_MRLR = matrix(0,n.rep,length(n.snp.vec))
beta_est_MRLR_inner = matrix(0,n.rep,length(n.snp.vec))
set.seed(i1)
for(m in 1:length(n.snp.vec)){
  #load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata")
  load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata")
  n.sub <- nrow(genotype_s)
  n.snp = n.snp.vec[m]
  genotype_s = genotype_s[,1:n.snp]
  
  beta_M = 0.15
  sigma_G  = 0.4
  #sigma_m = 1
  var_U = 1
  U = rnorm(n.sub,sd = sqrt(var_U))
  alpha_U <- 0.1
  beta_U <- 0.1
  sigma_u = alpha_U^2*var_U
  #sigma_e = 1-sigma_G-sigma_u
  sigma_e = 1-sigma_G
  #generate M phenotypes
  
  alpha_G = rnorm(n.snp,sd=sqrt(sigma_G/n.snp))
  #alpha_G[1] =0.1
  
  M_mat <- matrix(0,n.sub,n.rep)
  #G_value = genotype_s%*%alpha_G + alpha_U*U
  G_value = genotype_s%*%alpha_G 
  for(j in 1:n.rep){
    print(j)
    M_mat[,j] = G_value+rnorm(n.sub,sd = sqrt(sigma_e))
    #M_mat[,j] = G_value
    #+rnorm(n.sub,sd = sqrt(sigma_e))
  }
  # save(M_mat,file = "/data/zhangh24/MR_MA/result/simulation/prs/M_mat.rdata")
  
  
  #generate Y phenotypes
  #load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_Y_500k.rdata")
  load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_Y.rdata")
  genotype_s2 = genotype_s2[,1:n.snp]
  n.sub <- nrow(genotype_s2)
  U2 = rnorm(n.sub,sd = sqrt(var_U))
  n.rep = 2
  M_mat_inner = matrix(0,n.sub,n.rep)
  Y_mat <- matrix(0,n.sub,n.rep)
  G_value2 =  genotype_s2%*%alpha_G
  #+ alpha_U*U2
  #G_value2 =  genotype_s2%*%alpha_G 
  #+ alpha_U*U2
  #sigma_ey = sigma_G*beta_M^2/0.2-beta_M^2-beta_U^2*var_U
  sigma_ey = 0.2
  for(j in 1:n.rep){
    print(j)
    M =G_value2 + rnorm(n.sub,sd = sqrt(sigma_e))  
    M_mat_inner[,j] = M
    #M =G_value2 
    #+ rnorm(n.sub,sd = sqrt(sigma_e))  
    Y_mat[,j] = beta_M*M+rnorm(n.sub,sd=sqrt(sigma_ey))  
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
    model = lm(prs_y_mat~prs_m_mat)
    F = crossprod(G_value2)/n.snp/sigma_e
    
    beta_est_result[l,m]= coefficients(summary(model))[2,1]*(F+1)/F
    model  = lm(prs_y_mat~prs_m_mat_inner)
    beta_est_result_inner[l,m] = coefficients(summary(model))[2,1]
    
    
    #IVW method
    idx <- which(alpha_p<=5E-08)
    
    beta_est_IVW[l,m] = IVW_c(Gamma_est[idx],Gamma_sd[idx]^2,alpha_est[idx],alpha_sd[idx]^2)[1]
    beta_est_MRLR[l,m] = MRLR(Gamma_est[idx],Gamma_sd[idx]^2,alpha_est[idx],alpha_sd[idx]^2)[1]
    
    idx <- which(alpha_p_inner<=5E-08)
    beta_est_IVW_inner[l,m] = IVW_c_inner(Gamma_est[idx],Gamma_sd[idx]^2,alpha_est_inner[idx],alpha_sd_inner[idx]^2)[1]
    beta_est_MRLR_inner[l,m] = MRLR(Gamma_est[idx],Gamma_sd[idx]^2,alpha_est_inner[idx],alpha_sd_inner[idx]^2)[1]
    
  }
  
}
result= list(beta_est_result,
             beta_est_IVW,
             beta_est_MRLR,
             beta_est_result_inner,
             beta_est_IVW_inner,
             beta_est_MRLR_inner)
save(result,file  = paste0("/data/zhangh24/MR_MA/result/simulation/prs/beta_test_result_test_",i1))
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
