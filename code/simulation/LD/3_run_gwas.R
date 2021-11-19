args = commandArgs(trailingOnly = T)
#i represent the beta vector
#l represent the causal SNPs proportion
#v represent the pleotropic proportion
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
v = as.numeric(args[[3]])
#i1 represent the sub id, split 100 replciates into 10
i1 = as.numeric(args[[4]])
j = 22
num = 10
library(dplyr)
library(bc2)
library(data.table)
n.rep = n_rep = 100
start.end = startend(n.rep,num,i1)
start = start.end[1]
end = start.end[2]
#load the phenotpypes data and use plink to run
  cur.dir <- "/data/zhangh24/MR_MA/result/LD/"

  sid<-Sys.getenv('SLURM_JOB_ID')
  dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
  temp.dir = paste0('/lscratch/',sid,'/test/')
  system(paste0("cp ", cur.dir,"chr",j,".hm3.bed /lscratch/",sid,"/test/chr",j,".hm3.bed"))
  system(paste0("cp ", cur.dir,"chr",j,".hm3.bim /lscratch/",sid,"/test/chr",j,".hm3.bim"))
  system(paste0("cp ", cur.dir,"chr",j,".hm3.fam /lscratch/",sid,"/test/chr",j,".hm3.fam"))
  
  pheno_y = as.data.frame(fread(paste0(cur.dir,"y_pheno_plink_beta_",i,"_rho_",l,"_ple_",v,".phen")))
  pheno_y_sub = pheno_y[,c(1:2,(2+start):(2+end))]
  write.table(pheno_y_sub,file = paste0(temp.dir,"pheno_y_sub.phen"),row.names = F,col.names =F,quote=F )
  res <- system(paste0("/data/zhangh24/software/plink2_alpha --threads 2 --bfile /lscratch/",sid,"/test/chr",j,".hm3 --out ",temp.dir,"y_summary.out --glm omit-ref --pheno ",temp.dir,"pheno_y_sub.phen"))
  if(res==2){
    stop()
  }
  i_rep = 1
  sum.data = fread(paste0(temp.dir,"y_summary.out.PHENO",i_rep,".glm.linear"))
  sum.data.infor = sum.data[,1:6]
  sum.data.list = list()
  for(i_rep in 1:(end-start+1)){
    sum.data = fread(paste0(temp.dir,"y_summary.out.PHENO",i_rep,".glm.linear"))
    effect = sum.data[,c(9,10,12)]
    sum.data.list[[i_rep]] = effect
  }
  effect = bind_cols(sum.data.list)
  sum.data = cbind(sum.data.infor,effect)
  write.table(sum.data,file = paste0(cur.dir,"y_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_sub_",i1),row.names = F,col.names = T,quote=F)
  
  #system(paste0('rm -rf /lscratch/',sid,'/test/y_summary_chr_*'))
  
  pheno_m= as.data.frame(fread(paste0(cur.dir,"m_pheno_plink_beta_",i,"_rho_",l,"_ple_",v,".phen")))
  pheno_m_sub = pheno_m[,c(1:2,(2+start):(2+end))]
  write.table(pheno_m_sub,file = paste0(temp.dir,"pheno_m_sub.phen"),row.names = F,col.names =F,quote=F )
  # G_temp  = G[,idx[1]]
  # summary(lm(M_mat[,1]~G_temp))
  # 
  # head(pheno_m_sub)
  # summary(lm(pheno_m_sub[,3]~G_temp))
  # 
  #the SNP ID doens't represent the ref and alternative ID, that's the original 1KG ID
  #the reference and alternative ID should be based on the column
  res <- system(paste0("/data/zhangh24/software/plink2_alpha --threads 2 --bfile /lscratch/",sid,"/test/chr",j,".hm3 --out ",temp.dir,"m_summary.out --glm omit-ref --pheno ",temp.dir,"pheno_m_sub.phen"))
  if(res==2){
    stop()
  }
  i_rep = 1
  sum.data = fread(paste0(temp.dir,"m_summary.out.PHENO",i_rep,".glm.linear"))
  sum.data.infor = sum.data[,1:6]
  sum.data.list = list()
  for(i_rep in 1:(end-start+1)){
    sum.data = fread(paste0(temp.dir,"m_summary.out.PHENO",i_rep,".glm.linear"))
    effect = sum.data[,c(9,10,12)]
    sum.data.list[[i_rep]] = effect
  }
  effect = bind_cols(sum.data.list)
  sum.data = cbind(sum.data.infor,effect)
  write.table(sum.data,file = paste0(cur.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_sub_",i1),row.names = F,col.names = T,quote=F)
  
  
  

# for(i in 1:5){
#   fam <- data.frame(fread(paste0(cur.dir,eth[i],"/all_chr.tag.fam")))
#     for(j in 1:22){
#     write.table(fam, file = paste0(cur.dir,eth[i],"/chr",j,".tag.fam"),row.names = F,col.names = F,quote=F)
#   }
# }
# 
# 
