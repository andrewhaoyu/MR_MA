args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
i2 = as.numeric(args[[2]])
library(RcppArmadillo)
#library(withr)
#library(devtools)
#with_libpaths(new = "/home/zhangh24/R/x86_64-pc-linux-gnu-library/4.0/", install_github('andrewhaoyu/bc2'))
library(bc2, lib.loc ="/home/zhangh24/R/x86_64-pc-linux-gnu-library/4.0/")
FitLinearmodel <- function(y,x){
  model <- fastLm(X=x,y=y)
  if(is.na(coef(model)[2])){
    result <- c(0,1,1)
  }else{
    result <- coef(summary(model))[2,c(1,2,4)]  
  }
  return(result)
}
#run regression to obtain summary level statistics



if(i1==1){
  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata"))
  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_Y.rdata"))
  load("/data/zhangh24/MR_MA/result/simulation/prs/M_mat.rdata")
  load("/data/zhangh24/MR_MA/result/simulation/prs/Y_mat.rdata")
  
}else{
  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/noncau_genotype_M_",i1,".rdata"))
  load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/noncau_genotype_Y_",i1,".rdata"))
  load("/data/zhangh24/MR_MA/result/simulation/prs/M_mat.rdata")
  load("/data/zhangh24/MR_MA/result/simulation/prs/Y_mat.rdata")
}

n.train = 100000
genotype_m_train = genotype_s[1:n.train,]
M_mat_train = M_mat[1:n.train,]
n.rep = 1000
n.snp <- ncol(genotype_m_train)
size = 30
temp = startend(n.snp,size,i2)
start = temp[1]
end = temp[2]
n.snp.sub = end-start+1
alpha_G_est = matrix(0,n.snp.sub,n.rep)
sd_alpha_G = matrix(0,n.snp.sub,n.rep)
p_alpha_G = matrix(0,n.snp.sub,n.rep)
temp =1 
for(l in start:end){
  G_temp = cbind(1,genotype_m_train[,l])
  for(k in 1:n.rep){
    print(k)
    result_temp = FitLinearmodel(M_mat_train[,k],G_temp)
    alpha_G_est[temp,k]= result_temp[1]
    sd_alpha_G[temp,k]= result_temp[2]
    p_alpha_G[temp,k] = result_temp[3] 
  }
  temp = temp+1
}  


n.train = 100000
genotype_y_train = genotype_s2
Y_mat_train = Y_mat
n.rep = 1000
n.snp <- ncol(genotype_y_train)
size = 30
temp = startend(n.snp,size,i2)
start = temp[1]
end = temp[2]
n.snp.sub = end-start+1
gamma_G_est = matrix(0,n.snp.sub,n.rep)
sd_gamma_G = matrix(0,n.snp.sub,n.rep)
p_gamma_G = matrix(0,n.snp.sub,n.rep)
temp =1 
for(l in start:end){
  G_temp = cbind(1,genotype_y_train[,l])
  for(k in 1:n.rep){
    print(k)
    result_temp = FitLinearmodel(Y_mat_train[,k],G_temp)
    gamma_G_est[temp,k]= result_temp[1]
    sd_gamma_G[temp,k]= result_temp[2]
    p_gamma_G[temp,k] = result_temp[3] 
  }
  temp = temp+1
}  

result = list(alpha_G_est,sd_alpha_G,p_alpha_G,
              gamma_G_est,sd_gamma_G,p_gamma_G)

save(result,file = paste0("/data/zhangh24/MR_MA/result/simulation/prs/summary_gwas_",i1,"_",i2,".rdata"))
