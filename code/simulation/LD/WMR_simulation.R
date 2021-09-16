#Goal: Simulate data with realistic LD information
#Y = M\beta + G\theta + U\beta_u + error_y
#M = G\alpha + U\alpha_U + error_m
#h2_y  = var(\beta G\alpha+G\theta) = 0.4
#h2_m = var(\alphaG) = 0.4
#var(error_m+U\alpha_U) = 0.6
#var(error_y+U\beta_U) = 0.6
#causal SNPs proportion for M: 0.1, 0.01
#overlapping between pleotripic and non pleotropic 1, 0.5, 0.75
#i1 for beta
#i2 for ple
#i3 for rep
args = commandArgs(trailingOnly = T)
i1 = args[[1]]
i2 = args[[2]]
i3 = args[[3]]


setwd("/data/zhangh24/MR_MA/")
beta_vec = c(sqrt(0.4),0.3,0)
pleo_vec  = c(1,0.5,0.25)
n.snp = 500
beta = beta_vec[i1]
N = 60000
cau.pro = 0.2
n.cau = as.integer(n.snp*cau.pro)
h2_m = 0.4
h2_y = 0.4
sigma_alpha = h2_m/n.cau
sigma_theta = (h2_y-beta^2)/n.cau
alpha_u = 0.547
sigma_error_m = 1-h2_m-alpha_u^2
beta_u = 0.547
sigma_error_y = 1-h2_y-beta_u^2

set.seed(123)
idx.cau_m = sample(c(1:n.snp),n.cau)
#plotropic settings
pleosnp.pro = pleo_vec[i2]
n.cau.overlap = pleosnp.pro*n.cau
n.cau.specific = n.cau - n.cau.overlap
#pleotrpic snps proportion the same as causal snps
idx.cau_pleo = c(sample(idx.cau_m,n.cau.overlap),
                 sample(setdiff(c(1:n.cau),idx.cau_m),n.cau-n.cau.overlap))

alpha_G = rnorm(n.cau,mean = 0,sd = sqrt(sigma_alpha))
theta_G = rnorm(n.cau,mean = 0, sd = sqrt(sigma_theta))
ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
library(mr.raps)
library(Rfast)
library(MASS)
library(MESS)
R =ar1_cor(n.snp,0.5)
G1 = rmvnorm(N,mu = rep(0,n.snp),R)
G2 = rmvnorm(N,mu = rep(0,n.snp),R)
U1 = rnorm(N)
U2 = rnorm(N)
G1.cau = G1[,idx.cau_m]
G2.cau = G2[,idx.cau_m]
G1.pleo = G1[,idx.cau_pleo]
G2.pleo = G2[,idx.cau_pleo]
ldscore = rep(sum(R[2:15,]^2),n.snp)

#G1 to obtain sum data for Y
#G2 to obtain sum data for M
n.rep = 25
beta_est = rep(0,n.rep)
beta_cover = rep(0,n.rep)
beta_se = rep(0,n.rep)


beta_est_Raps = rep(0,n.rep)
beta_cover_Raps = rep(0,n.rep)
beta_se_Raps = rep(0,n.rep)
beta_est_IVW = rep(0,n.rep)
beta_cover_IVW = rep(0,n.rep)
beta_se_IVW = rep(0,n.rep)
beta_est_egger = rep(0,n.rep)
beta_cover_egger = rep(0,n.rep)
beta_se_egger  = rep(0,n.rep)
beta_est_median = rep(0,n.rep)
beta_cover_median = rep(0,n.rep)
beta_se_median = rep(0,n.rep)
library(MendelianRandomization)
library(susieR)
for(k in 1:n.rep){
  print(k)
  error_m = rnorm(N,sd = sqrt(sigma_error_m))
  error_y = rnorm(N,sd = sqrt(sigma_error_y))
  M1 = G1.cau%*%alpha_G+U1*alpha_u+error_m
  Y1 = M1%*%beta + G1.pleo%*%theta_G+U1*beta_u + error_y
  #Y1 = M1%*%beta +U1*beta_u + error_y
  error_m = rnorm(N,sd = sqrt(sigma_error_m))
  M2 = G2.cau%*%alpha_G+U2*alpha_u+error_m
  
  
  sumstats <- univariate_regression(G1, Y1)
  
  Gamma = sumstats$betahat
  se_Gamma = sumstats$sebetahat
  
  sumstats <- univariate_regression(G2, M2)
  alpha = sumstats$betahat
  se_alpha = sumstats$sebetahat
  p_alpha = 2*pnorm(-abs(alpha/se_alpha),lower.tail = T)
  #p_Gamma = 2*pnorm(-abs(Gamma/sqrt(var_Gamma)),lower.tail = T)
  Myclumping <- function(R,p){
    n.snp = ncol(R)
    
    
    #keep snps for clumpinp
    keep.ind = c(1:n.snp)
    #remove snps due to clumping
    remove.ind  = NULL
    #select snp ind
    select.ind = NULL
    temp = 1
    while(length(keep.ind)>0){
      # print(temp)
      p.temp = p[keep.ind]
      #select top snp
      top.ind = which.min(p.temp)
      select.ind = c(select.ind,keep.ind[top.ind])
      #print(keep.ind[top.ind])
      #tempory correlation
      R.temp = R[keep.ind[top.ind],]
      idx.remove = which(R.temp>=0.01)
      #take out
      remove.ind= c(remove.ind,idx.remove)
      keep.ind = setdiff(keep.ind,remove.ind)
      temp = temp+1
    }
    result = data.frame(select.ind,p[select.ind])
    return(result)
  }
  clump.snp = Myclumping(R,p_alpha)
  select.id = clump.snp[clump.snp$p.select.ind.<=5E-08,1]
  alpha_select =alpha[select.id]
  se_alpha_select = se_alpha[select.id]
  Gamma_select = Gamma[select.id]
  se_Gamma_select = se_Gamma[select.id]
  
  MRInputObject <- mr_input(bx = alpha_select,
                            bxse = se_alpha_select,
                            by = Gamma_select,
                            byse = se_Gamma_select)
  IVWObject <- mr_ivw(MRInputObject,
                      model = "default",
                      robust = FALSE,
                      penalized = FALSE,
                      correl = FALSE,
                      weights = "simple",
                      psi = 0,
                      distribution =
                      "normal",
                      alpha = 0.05)
  
  beta_est_IVW[k] = IVWObject$Estimate
  beta_cover_IVW[k] = ifelse(IVWObject$CILower<=beta&
                               IVWObject$CIUpper>=beta,1,0)
  beta_se_IVW[k] = IVWObject$StdError
  
  
  EggerObject <- mr_egger(
    MRInputObject,
    robust = FALSE,
    penalized = FALSE,
    correl = FALSE,
    distribution = "normal",
    alpha = 0.05
  )
    
  beta_est_egger[k] = EggerObject$Estimate
  beta_cover_egger[k] = ifelse(EggerObject$CILower.Est<=beta&
                                 EggerObject$CIUpper.Est>=beta,1,0)
  beta_se_egger[k] = EggerObject$StdError.Est
  
  MedianObject <- mr_median(
    MRInputObject,
    weighting = "weighted",
    distribution = "normal",
    alpha = 0.05,
    iterations = 10000,
    seed = 314159265
  )
  
  
  beta_est_median[k] = MedianObject$Estimate
  beta_cover_median[k] = ifelse(MedianObject$CILower<=beta&
                               MedianObject$CIUpper>=beta,1,0)
  beta_se_median[k] = MedianObject$StdError
  
  
  
  
  raps_result <- mr.raps(data = data.frame(beta.exposure = alpha_select,
                                           beta.outcome = Gamma_select,
                                           se.exposure = se_alpha_select,
                                           se.outcome = se_Gamma_select),
                         diagnostics = F)
  beta_est_Raps[k] = raps_result$beta.hat
  beta_cover_Raps[k] = ifelse(raps_result$beta.hat-1.96*raps_result$beta.se<=beta&
                                raps_result$beta.hat+1.96*raps_result$beta.se>=beta,1,0)
  beta_se_Raps[k] = raps_result$beta.se
 # se_Gamma = sqrt(var_Gamma)
  #se_alpha = sqrt(var_alpha)
  
   #R.select = R[select.id,select.id]
  MR_result <- WMRFun(Gamma,se_Gamma,
                                 alpha,se_alpha,
                                 ldscore,R)
    # MRWeight(Gamma = sumGamma,
    #                     var_Gamma = var_Gamma,
    #                     alpha = sumalpha,
    #                     var_alpha = var_alpha,
    #                     R = R)
  beta_est[k] = MR_result[1]
  beta_cover[k] = ifelse(MR_result[3]<=beta&MR_result[4]>=beta,1,0)
  beta_se[k] = MR_result[2]
  
}

mean.result = data.frame(
  beta_est,beta_est_IVW,beta_est_egger,beta_est_median,beta_est_Raps
)
colnames(mean.result) = c("WMR","IVW","MR-Egger","MR-median","MRRAPs")

se.result = data.frame(
  beta_se,beta_se_IVW,beta_se_egger,beta_se_median,beta_se_Raps
)
colnames(se.result) = c("WMR","IVW","MR-Egger","MR-median","MRRAPs")

cover.result = data.frame(
  beta_cover,beta_cover_IVW,beta_cover_egger,beta_cover_median,beta_cover_Raps
)
colnames(cover.result) = c("WMR","IVW","MR-Egger","MR-median","MRRAPs")

result = list(mean.result,se.result,cover.result)

save(result,file = paste0("./result/simulation/LD_simulation_test/result_",i1,"_",i2,"_",i3,".rdata"))
