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
sigma_e = 1-sigma_G-sigma_u
#generate M phenotypes
alpha_G = rnorm(n.snp,sd=sqrt(sigma_G/n.snp))
n.rep = 1000
M_mat <- matrix(0,n.sub,n.rep)
G_value = genotype_s%*%alpha_G + alpha_U*U
for(j in 1:n.rep){
  print(j)
  M_mat[,j] = G_value+rnorm(n.sub,sd = sqrt(sigma_e))
    }
save(M_mat,file = "/data/zhangh24/MR_MA/result/simulation/prs/M_mat.rdata")



#generate Y phenotypes
n.sub <- nrow(genotype_s2)
U2 = rnorm(n.sub,sd = sqrt(var_U))
n.rep = 1000
Y_mat <- matrix(0,n.sub,n.rep)
G_value2 =  genotype_s2%*%alpha_G + alpha_U*U2
sigma_ey = 1-beta_M^2-beta_U^2*var_U
for(j in 1:n.rep){
  print(j)
  M =G_value2 +rnorm(n.sub,sd = sqrt(sigma_e))  
  Y_mat[,j] = beta_M*M+beta_U*U2+rnorm(n.sub,sd=sqrt(sigma_ey))  
}

save(Y_mat,file = "/data/zhangh24/MR_MA/result/simulation/prs/Y_mat.rdata")

