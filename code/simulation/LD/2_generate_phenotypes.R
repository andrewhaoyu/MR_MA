
library(bigsnpr)
library(dplyr)
library(MASS)
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
cau_vec = c(0.05,0.01,0.001)
beta_vec = c(0,0.2)

pleo_vec= c(0,0.3,0.8)



snp.infor.update = left_join(map,snp.infor.subset,by="SNP")
MAF = snp.infor.update[,"MAF",drop=F]
#hapmap3 contains 1217311 SNPs
#assume total h2 as 0.4
#CHR 22 has probablity 17268/1217311*0.4= 0.0057
h2_y = h2_m = 0.057
set.seed(777)
n.snp = nrow(map)
n.rep = 100
cau_snp_list = list()
sample_size_vec = c(60000)
  for(l in 1:3){
    cau = cau_vec[l]
    n.cau = as.integer(n.snp*cau)
    idx.cau_m <- sample(c(1:n.snp),n.cau)
    sigma_alpha = h2_m/n.cau
    
    alpha_G = rnorm(n.cau,mean = 0,sd = sqrt(sigma_alpha))
    
    rho = 0.5
    sigma_error_m  = 1-h2_m
    sigma_error_y = 1-h2_y
    cov_my = sqrt(sigma_error_m*sigma_error_y)*rho
    Sigma = matrix(c(sigma_error_m,cov_my,cov_my,sigma_error_y),2,2)
    n.sub = nrow(G)
    
    alpha = rep(0,n.snp)
    alpha[idx.cau_m] = alpha_G/sqrt(2*MAF[idx.cau_m,1]*(1-MAF[idx.cau_m,1]))
    Galpha = big_prodVec(G, alpha, ind.row = rows_along(G), ind.col = cols_along(G))
    
    
    cau_snp_list[[l]] = data.frame(snp.infor.update[idx.cau_m,],
                                   alpha = alpha[idx.cau_m])
    
    for(v in 1:3){
      pleo.pro = pleo_vec[v]
      n.cau.overlap = as.integer(pleo.pro*n.cau)
      n.cau.specific = n.cau - n.cau.overlap
      idx.cau_y = c(sample(idx.cau_m,n.cau.overlap),
                    sample(setdiff(c(1:n.snp),idx.cau_m),n.cau-n.cau.overlap))
      
      for(i in 1:length(beta_vec)){
      beta = beta_vec[i]
      sigma_theta = (h2_y-beta^2*h2_m)/n.cau
      
      theta_G = rnorm(n.cau,mean = 0, sd = sqrt(sigma_theta))
      # cau_snp_y_list[[v]] = data.frame(snp.infor.update[idx.cau_y,],
      #                                  alpha = alpha[idx.cau_m])
      #if(n.cau.overlap!=0){
      #sample causal SNPs for y by controlling pleotropic
      theta = rep(0,n.snp)
      theta[idx.cau_y] = theta_G/sqrt(2*MAF[idx.cau_y,1]*(1-MAF[idx.cau_y,1]))
      Gtheta = big_prodVec(G,theta,ind.row = rows_along(G),ind.col = cols_along(G))
      #pleotrpic snps proportion the same as causal snps
      M_mat = matrix(0,n.sub,n.rep)
      Y_mat = matrix(0,n.sub,n.rep)
      for(i_rep in 1:n.rep){
        error = mvrnorm(n.sub,mu = c(0,0),Sigma = Sigma)  
        M_mat[,i_rep] = Galpha+error[,1]
        Y_mat[,i_rep] = beta*M_mat[,i_rep]+Gtheta+error[,2]
      }
      # G_temp  = G[,idx[1]]
      # summary(lm(M_mat[,1]~G_temp))
      # 
      #create phenotypes files for plink
      pheno_M <- fam[,1:2]
      pheno_M2 <- fam[,1:2]
      pheno_Y = fam[,1:2]
      for(m in 1:length(sample_size_vec)){
        #first set for selection
        temp_mat = M_mat
        temp_mat[40001:120000,] = NA
        pheno_M = cbind(pheno_M,temp_mat)
        #second set for exposure
        temp_mat = M_mat
        temp_mat[c(1:40000,80000:120000),] = NA
        pheno_M2 = cbind(pheno_M2,temp_mat)
        
        #third set for outcome
        Y_mat[1:80000,] = NA
        pheno_Y = cbind(pheno_Y,Y_mat)
      }
      write.table(pheno_Y,file = paste0(cur.dir,"y_pheno_plink_beta_",i,"_rho_",l,"_ple_",v,".phen"),row.names = F,col.names = F,quote=F)
      write.table(pheno_M,file = paste0(cur.dir,"m_pheno_plink_beta_",i,"_rho_",l,"_ple_",v,".phen"),row.names = F,col.names = F,quote=F)
      write.table(pheno_M2,file = paste0(cur.dir,"m2_pheno_plink_beta_",i,"_rho_",l,"_ple_",v,".phen"),row.names = F,col.names = F,quote=F)
      
    }

    
  
  }
  
}
save(cau_snp_list,file = paste0(cur.dir,"cau_snp_list.rdata"))  

