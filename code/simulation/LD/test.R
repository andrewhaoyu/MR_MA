set.seed(i1)
library(susieR)
setwd("/data/zhangh24/MR_MA/")
source("./code/simulation/functions/MR_weight_new.R")
setwd("/Users/zhangh24/GoogleDrive/MR_MA/result/simulation/LD")
load("example_G.rdata")
library(MASS)
n.sub = 60000
G1 = G.temp[1:n.sub,]
G2 = G.temp[(n.sub+1):(2*n.sub),]
n.snp = ncol(G1)
b <- rep(0,1000)
set.seed(666)
n.cau = 6
idx.cau = c(250,245,500,505,700,710)
R[idx.cau,idx.cau]
b[idx.cau] = rnorm(n.cau,0,sd = 0.04)
idx <- which(b!=0)
plot(b, pch=16, ylab='effect size')
MAF = colSums(G)/(nrow(G))

beta_M = 0.15
sigma_m = 0.6
beta_M = 0.15
rho = 0.3
sigma_y = 1-beta_M^2
sigma_ym = sigma_y*sigma_m*rho
Sigma = matrix(c(sigma_y,sigma_ym,sigma_ym,sigma_m),2,2)
error = mvrnorm(n.sub,mu = c(0,0),Sigma = Sigma) 
M1 = G1%*%b+error[,1]
Y1 = M1*beta_M+error[,2]
M2 = G2%*%b+rnorm(n.sub,mean = 0, sd = sqrt(sigma_m))
sumstats_M <- univariate_regression(G2, M2)
sumstats_Y <- univariate_regression(G1, Y1)
z_scores <- sumstats_M$betahat / sumstats_M$sebetahat
p.value = 2*pnorm(-abs(z_scores),lower.tail=T)

cau = rep("Not causal SNPs",n.snp)
cau[idx.cau] = "Causal SNPs"
SNP.index = c(1:n.snp)
plot.data = data.frame(SNP.index,p = -log10(p.value),cau)
library(ggplot2)

ggplot(plot.data,aes(x= SNP.index,
                     y = p,col=cau))+
  geom_point()+
  scale_colour_manual(values =c("red","black"))+
  ylab("-log10(p)")+
  xlab("SNP")+
  theme(legend.title = element_blank())+
  theme_Publication()
library(ggplot2)

#sample subjects to calculate correlations
idx <- sample(c(1:nrow(G)),3000)
G_sub = G.temp[idx,]
R = cor(G_sub)
R2 = R^2
R[idx.cau,idx.cau]
R_temp = R
p_temp = p.value

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




clump = Myclumping(R2,p.value)

select.ind = clump[,1]


idx = which(clump[,2]<=5E-08)
select.ind = clump[idx,1]

Gamma = sumstats_Y$betahat[select.ind]
var_Gamma = sumstats_Y$sebetahat[select.ind]^2
gamma = sumstats_M$betahat[select.ind]
var_gamma = sumstats_Y$sebetahat[select.ind]^2
library(mr.raps)
raps_result <- mr.raps(data = data.frame(beta.exposure = gamma,
                                         beta.outcome = Gamma,
                                         se.exposure = sqrt(var_gamma),
                                         se.outcome = sqrt(var_Gamma)),
                       diagnostics = F)
raps_result$beta.hat
raps_result$beta.se





IVW_c_result<- IVW_c(Gamma,var_Gamma,
                    gamma,var_gamma)


fitted_rss <- susie_rss(z_scores, R, L = 10,
                        estimate_residual_variance = TRUE, 
                        estimate_prior_variance = TRUE)

susie_plot(fitted_rss, y="PIP", b=b)
library(dplyr)

FindIV <- function(fitted_rss){
  fitted_rss_temp = as.data.frame(summary(fitted_rss)$vars) %>% 
    filter(cs!=-1) %>% 
    group_by(cs) %>% 
    filter(variable_prob == max(variable_prob)) %>% 
    select(variable)
  select = as.vector(fitted_rss_temp$variable)
  
  return(select)
}



select.ind = FindIV(fitted_rss)

Gamma = sumstats_Y$betahat[select.ind]
var_Gamma = sumstats_Y$sebetahat[select.ind]^2
gamma = sumstats_M$betahat[select.ind]
var_gamma = sumstats_Y$sebetahat[select.ind]^2

R.select= R[select.ind,select.ind]
R_keep = R


R = R.select
MR_weight_result = MRWeight(Gamma,var_Gamma,
                gamma,var_gamma,R.select)



result = list(raps_result,IVW_c_result,MR_weight_result)