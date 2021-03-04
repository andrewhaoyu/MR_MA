
library(bigsnpr)
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
j = 22
#snp_readBed(paste0(cur.dir,"chr",j,".hm3.bed"))
obj.bigSNP <- snp_attach(paste0(cur.dir,"chr",j,".hm3.rds"))
G  = obj.bigSNP$genotypes
fam  = obj.bigSNP$fam
map = obj.bigSNP$map
setwd("/data/zhangh24/multi_ethnic/")
load("/data/zhangh24/MR_MA/result/LD/chr22_snp_infor.rdata")
colnames(map)[2] <- "SNP"
cau_vec = c(0.01,0.001,5E-04)
snp.infor.update = left_join(map,snp.infor.subset,by="SNP")
MAF = snp.infor.update %>% select(EUR)
#hapmap3 contains 1217311 SNPs
#assume total h2 as 0.4
#CHR 22 has probablity 17268/1217311*0.4= 0.0057
h2 = 0.0057
set.seed(666)
n.snp = nrow(map)
n.rep = 100
sigma_m = 1-h2
beta_M = 0.15
rho = 0.3
sigma_y = 1-beta_M^2
sigma_ym = sigma_y*sigma_m*rho
Sigma = matrix(c(sigma_y,sigma_ym,sigma_ym,sigma_m),2,2)
n.sub = nrow(G)
library(MASS)
sample_size_vec = c(60000)
for(l in 1:3){
  cau = cau_vec[l]
  n.cau = as.integer(n.snp*cau)
  u = rnorm(n.cau,mean = 0,sd  = sqrt(h2/n.cau))
  idx <- sample(c(1:n.snp),n.cau)
  alpha = rep(0,n.snp)
  alpha[idx] = u/sqrt(2*MAF[idx,1]*(1-MAF[idx,1]))
  Galpha = big_prodVec(G, alpha, ind.row = rows_along(G), ind.col = cols_along(G))
  M_mat = matrix(0,n.sub,n.rep)
  Y_mat = matrix(0,n.sub,n.rep)
  for(i_rep in 1:n.rep){
    error = mvrnorm(n.sub,mu = c(0,0),Sigma = Sigma)  
    M_mat[,i_rep] = Galpha+error[,1]
    Y_mat[,i_rep] = beta_M*Galpha+error[,2]
  }
  #create phenotypes files for plink
  pheno_M <- fam[,1:2]
  pheno_Y = fam[,1:2]
  for(m in 1:length(sample_size_vec)){
    M_mat[1:sample_size_vec[m],] = NA
    pheno_M = cbind(pheno_M,M_mat)
    Y_mat[(60000+1):(60000+sample_size_vec[m]),] = NA
    pheno_Y = cbind(pheno_Y,Y_mat)
  }
  write.table(pheno_M,file = paste0(cur.dir,"y_pheno_plink_rho_",l,".phen"),row.names = F,col.names = F,quote=F)
  write.table(pheno_Y,file = paste0(cur.dir,"m_pheno_plink_rho_",l,".phen"),row.names = F,col.names = F,quote=F)
}
