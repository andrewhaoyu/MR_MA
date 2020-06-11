#generate phenotypes for the two-sample MR analysis

set.seed(666)
load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata")
n.sub <- nrow(genotype_s)
n.snp = 5000
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
n.rep = 100
M_mat <- matrix(0,n.sub,n.rep)
#G_value = genotype_s%*%alpha_G + alpha_U*U
G_value = genotype_s%*%alpha_G 
for(j in 1:n.rep){
  print(j)
 # M_mat[,j] = G_value+rnorm(n.sub,sd = sqrt(sigma_e))
  M_mat[,j] = G_value
  #+rnorm(n.sub,sd = sqrt(sigma_e))
    }
save(M_mat,file = "/data/zhangh24/MR_MA/result/simulation/prs/M_mat.rdata")



#generate Y phenotypes
load("/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_Y.rdata")
n.sub <- nrow(genotype_s2)
#U2 = rnorm(n.sub,sd = sqrt(var_U))
n.rep = 100
Y_mat <- matrix(0,n.sub,n.rep)
#G_value2 =  genotype_s2%*%alpha_G + alpha_U*U2
G_value2 =  genotype_s2%*%alpha_G 
#+ alpha_U*U2
#sigma_ey = sigma_G*beta_M^2/0.2-beta_M^2-beta_U^2*var_U
sigma_ey = 0.01
for(j in 1:n.rep){
  print(j)
  #M =G_value2 + rnorm(n.sub,sd = sqrt(sigma_e))  
  M =G_value2 
  #+ rnorm(n.sub,sd = sqrt(sigma_e))  
  #Y_mat[,j] = beta_M*M+rnorm(n.sub,sd=sqrt(sigma_ey))  
  Y_mat[,j] = beta_M*M
  #+rnorm(n.sub,sd=sqrt(sigma_ey))  
  #Y_mat[,j] = beta_M*M+beta_U*U2+rnorm(n.sub,sd=sqrt(sigma_ey))  
}

beta_est = rep(0,n.rep)
for(k in 1:n.rep){
  print(k)
  model1 <- lm(M_mat[,k]~genotype_s[,64])
  alpha_est = coef(summary(model1))[2,1]
  model2 <- lm(Y_mat[,k]~genotype_s2[,64])
  Gamma_est = coef(summary(model2))[2,1]

  beta_est[k] = Gamma_est/alpha_est
}




save(Y_mat,file = "/data/zhangh24/MR_MA/result/simulation/prs/Y_mat.rdata")

