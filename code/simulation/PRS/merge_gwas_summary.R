#merge the GWAS summary statistics from R
n.snp <- 150000
n.rep <- 1000

alpha_mat <- matrix(0,n.snp,n.rep)
alpha_sd_mat <- matrix(0,n.snp,n.rep)
alpha_p_mat <- matrix(0,n.snp,n.rep)
gamma_mat <- matrix(0,n.snp,n.rep)
gamma_sd_mat <- matrix(0,n.snp,n.rep)
gamma_p_mat <- matrix(0,n.snp,n.rep)

total <- 0
for(i1 in 1:30){
  for(i2 in 1:30){
   load(paste0("/data/zhangh24/MR_MA/result/simulation/prs/summary_gwas_",i1,"_",i2,".rdata")) 
    temp = nrow(result[[1]])
    alpha_mat[total+(1:temp),] = result[[1]]
    alpha_sd_mat[total+(1:temp),] = result[[2]]
    alpha_p_mat[total+(1:temp),] = result[[3]]
    gamma_mat[total+(1:temp),]  = result[[4]]
    gamma_sd_mat[total+(1:temp),] = result[[5]]
    gamma_p_mat[total+(1:temp),] = result[[6]]
    total = total+temp
  }
}

result <- list(alpha_mat,alpha_sd_mat,alpha_p_mat,
               gamma_mat,gamma_sd_mat,gamma_p_mat)
save(result,file = "/data/zhangh24/MR_MA/result/simulation/prs/summary_gwas_our.rdata")
