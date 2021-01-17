pthres <- c(5E-08,1E-07,5E-07,1E-06,5E-06,1E-05,5E-05,1E-04,1E-03,1E-02,1E-01)

load("/data/zhangh24/MR_MA/result/simulation/prs/summary_gwas_our.rdata")
alpha_est_mat = result[[1]]
alpha_sd_mat = result[[2]]
alpha_p_mat = result[[3]]
gamma_est_mat = result[[4]]
gamma_sd_mat = result[[5]]
gamma_p_mat = result[[6]]
load("/data/zhangh24/MR_MA/result/simulation/prs/M_mat.rdata")
n.train = 100000
M_mat = M_mat[1:n.train,]
# result_new <- list(alpha_est_mat[idx,1:10],gamma_est_mat[idx,1:10])
# save(result_new,file = "/data/zhangh24/MR_MA/result/simulation/prs/summary_gwas_06.rdata")

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
  
  return(c(ratio_ivw,ratio_ivw_var,cover,
           coef_low,coef_high,p_value))
}
Meta = function(coef_vec,var_vec){
  meta_var = (sum(1/var_vec))^-1
  meta_coef = meta_var*sum(coef_vec/var_vec)
  return(c(meta_coef,meta_var))
}

# n.rep = ncol(alpha_est_mat)
# IVW_best_est = rep(0,n.rep)
# IVW_best_low_est = rep(0,n.rep)
# IVW_best_high_est = rep(0,n.rep)
# IVW
# for(l in 1:n.rep){
#   #idx = 1:5000
#   Gamma = gamma_est_mat[idx,l]
#   var_Gamma = gamma_sd_mat[idx,l]^2
#   alpha = alpha_est_mat[idx,l]
#   var_alpha = alpha_sd_mat[idx,l]^2
# 
#   result_IVW=  IVW_c(Gamma,var_Gamma,
#                       alpha,var_alpha)
#   IVW_best_est[l] = result_IVW[1]
#   IVW_best_low_est[l] = result_IVW[4]
#   IVW_best_high_est[l] = result_IVW[5]
#    #IVW_cover[l] = result_IVW[3]
#   # 
# }





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
#   cover <- ifelse((beta_M>=coef_low&
#                      beta_M<=coef_high),1,0)
#   
#   return(list(coef_est,coef_low,coef_high,cover))
# }









load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/prs_result_combined_twostage.rdata"))
prs_m_mat <- prs_result[[1]]
prs_y_mat <- prs_result[[2]]
n.rep = 1000
#beta_M = mean(best_prs_est)
# MRLR_est = rep(0,n.rep)
# MRLR_low_est = rep(0,n.rep)
# MRLR_high_est = rep(0,n.rep)
# MRLR_cover = rep(0,n.rep)
# 
# MRLR_PRS_est = rep(0,n.rep)
# MRLR_PRS_low_est = rep(0,n.rep)
# MRLR_PRS_high_est = rep(0,n.rep)
# MRLR_PRS_cover = rep(0,n.rep)

# IVW_p = rep(0,n.rep)
# IVW_PRS_est = rep(0,n.rep)
# IVW_PRS_low_est = rep(0,n.rep)
# IVW_PRS_high_est = rep(0,n.rep)
# IVW_PRS_cover = rep(0,n.rep)
# IVW_PRS_p = rep(0,n.rep)
IVW_est = rep(0,n.rep)
IVW_low_est = rep(0,n.rep)
IVW_high_est = rep(0,n.rep)
IVW_p = rep(0,n.rep)
IVW_cover = rep(0,n.rep)
beta_M = 0.15

MRPRS_est = rep(0,n.rep)
MRPRS_low_est = rep(0,n.rep)
MRPRS_high_est = rep(0,n.rep)
MRPRS_cover = rep(0,n.rep)
MRPRS_p = rep(0,n.rep)

for(l in 1:n.rep){
  print(l)
  #IVW method restricted to genome-wide SNPs
  idx <- which(alpha_p_mat[,l]<=5E-8)
  Gamma = gamma_est_mat[idx,l]
  var_Gamma = gamma_sd_mat[idx,l]^2
  alpha = alpha_est_mat[idx,l]
  var_alpha = alpha_sd_mat[idx,l]^2
  result_IVW = IVW_c(Gamma,var_Gamma,
        alpha,var_alpha)
  IVW_est[l] = result_IVW[1]
  IVW_low_est[l] = result_IVW[4]
  IVW_high_est[l] = result_IVW[5]
  IVW_cover[l] = result_IVW[3]
  IVW_p[l] = result_IVW[6]
  #MR-PRS method restricted to best C+T SNPs
  model_m = lm(M_mat[,l]~prs_m_mat[,l])
  sigma_m = summary(model_m)$sigma
  idx = which(alpha_p_mat[,l]<=1E-03)
  Q = length(idx)
  F =   crossprod(prs_m_mat[,l])/sigma_m/Q
  model <- lm(prs_y_mat[,l]~prs_m_mat[,l])
    MRPRS_est[l] <- coefficients(summary(model))[2,1]*(F+1)/F
    
    temp_ci <- confint(model)
    MRPRS_low_est[l] <- temp_ci[2,1]*(F+1)/F
    MRPRS_high_est[l] <- temp_ci[2,2]*(F+1)/F
    MRPRS_cover[l] = ifelse((beta_M>=temp_ci[2,1])&
                              (beta_M<=temp_ci[2,2]),1,0)
    MRPRS_p[l] <- coefficients(summary(model))[2,4]
    
}
  # result_MRLR = MRLR(Gamma,var_Gamma,
  #                    alpha,var_alpha)
  # MRLR_est[l] = as.numeric(result_MRLR[1])
  # MRLR_low_est[l] = as.numeric(result_MRLR[2])
  # MRLR_high_est[l] = as.numeric(result_MRLR[3])
  # MRLR_cover[l] = as.numeric(result_MRLR[4])
  
  
  # idx <- which(alpha_p_mat[,l]<=1E-03)
  # Gamma = gamma_est_mat[idx,l]
  # var_Gamma = gamma_sd_mat[idx,l]^2
  # alpha = alpha_est_mat[idx,l]
  # var_alpha = alpha_sd_mat[idx,l]^2
  # result_IVW = IVW_c(Gamma,var_Gamma,
  #                    alpha,var_alpha)
  # IVW_PRS_est[l] = result_IVW[1]
  # IVW_PRS_low_est[l] = result_IVW[4]
  # IVW_PRS_high_est[l] = result_IVW[5]
  # IVW_PRS_cover[l] = result_IVW[3]
  # IVW_PRS_p[l] = result_IVW[6]
  # result_PRS_MRLR = MRLR(Gamma,var_Gamma,
  #                        alpha,var_alpha)
  # MRLR_PRS_est = as.numeric(result_PRS_MRLR[1])
  # MRLR_PRS_low_est = as.numeric(result_PRS_MRLR[2])
  # MRLR_PRS_high_est = as.numeric(result_PRS_MRLR[3])
  # MRLR_PRS_cover = as.numeric(result_PRS_MRLR[4])
  
  #PRS_est[l] = crossprod(Gamma,alpha)/crossprod(alpha,alpha)
  
}

est = c(mean(IVW_est),mean(MRLR_est),
        mean(IVW_PRS_est),mean(MRLR_PRS_est),
        mean(best_prs_est),
        mean(IVW_best_est),
        mean(prs_est))
est_low = c(mean(IVW_low_est),mean(MRLR_low_est),
            mean(IVW_PRS_low_est),mean(MRLR_PRS_low_est),
            mean(best_prs_low),
            mean(IVW_best_low_est),
            mean(prs_low_est))
est_high = c(mean(IVW_high_est),mean(MRLR_high_est),
            mean(IVW_PRS_high_est),mean(MRLR_PRS_high_est),
            mean(best_prs_high),
            mean(IVW_best_high_est),
            mean(prs_high_est))
method = c(
           )
data <- data.frame(method,est,est_low,est_high)
data$method <- as.factor(data$method)
levels(data$method) <- c("Two-stage regression (best PRS)",
                         "IVW true causal SNPS",
                         "IVW (P<5E-08)","MR-weighted (P<5E-08)",
                         "IVW (P<1E-03)"," MR-weighted (P<1E-03)",
                         "Two-stage regression (P<1E-03)")
save(data,file = "/data/zhangh24/MR_MA/result/simulation/prs/data_for_plot.rdata")
library(ggplot2)
 ggplot(data) +
  #theme_Publication()+
  geom_bar(aes(x=method, y=est), stat="identity", fill="royalblue", alpha=0.7) +
  geom_errorbar(aes(x=name, ymin=est_low, ymax=est_high, width=0.2), 
                colour="firebrick2", alpha=0.9, size=0.8) +
  coord_flip() +
  # ylim(-5, 15) +
  labs(x = 'Variables', y = expression("95% CI"), 
       title = 'MR method comparasion') +
  theme(plot.title=element_text(hjust=0.5, size=30),
        plot.subtitle=element_text(size=20),
        axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0), hjust=0.5, size=15),
        axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0), vjust=0.5, size=15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15))
