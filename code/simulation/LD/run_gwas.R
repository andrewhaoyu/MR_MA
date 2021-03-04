
l = as.numeric(args[[1]])
#i1 represent the sub id, split 100 replciates into 10
i1 = as.numeric(args[[2]])
j = 22
num = 10
library(dplyr)
library(bc2)
library(data.table)
n_rep = 100
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
  
  pheno_y = as.data.frame(fread(paste0(cur.dir,"y_pheno_plink_rho_",l,".phen")))
  
  pheno_y_sub = pheno_y[,c(1:2,(2+start):(2+end))]
  write.table(pheno_y_sub,file = paste0(temp.dir,"pheno_y_sub_",i1,".phen"),row.names = F,col.names =F,quote=F )
  res <- system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /lscratch/",sid,"/test/chr",j,".hm3 --out ",temp.dir,"y_summary_chr_",j,"_rho_",l,".out --linear --all-pheno --allow-no-sex --pheno ",temp.dir,"pheno_y_sub_",i1,".phen"))
  if(res==2){
    stop()
  }
  i_rep = 1
  sum.data = fread(paste0(temp.dir,"y_summary_chr_",j,"_rho_",l,".out.P",i_rep,".assoc.linear"))
  sum.data.infor = sum.data[,1:6]
  sum.data.list = list()
  for(i_rep in 1:(end-start+1)){
    sum.data = fread(paste0(temp.dir,"y_summary_chr_",j,"_rho_",l,".out.P",i_rep,".assoc.linear"))
    effect = sum.data[,c(7:9)]
    sum.data.list[[i_rep]] = effect
  }
  effect = bind_cols(sum.data.list)
  sum.data = cbind(sum.data.infor,effect)
  write.table(sum.data,file = paste0(cur.dir,"y_summary_chr_",j,"_rho_",l,"_sub_",i1),row.names = F,col.names = T,quote=F)
  
  #system(paste0('rm -rf /lscratch/',sid,'/test/y_summary_chr_*'))
  
  pheno_m= as.data.frame(fread(paste0(cur.dir,"m_pheno_plink_rho_",l,".phen")))
  pheno_m_sub = pheno_m[,c(1:2,(2+start):(2+end))]
  write.table(pheno_m_sub,file = paste0(temp.dir,"pheno_m_sub_",i1,".phen"),row.names = F,col.names =F,quote=F )
  
  res <- system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /lscratch/",sid,"/test/chr",j,".hm3 --out ",temp.dir,"m_summary_chr_",j,"_rho_",l,".out --linear --all-pheno --allow-no-sex --pheno ",cur.dir,"y_pheno_plink_rho_",l,".phen"))
  if(res==2){
    stop()
  }
  i_rep = 1
  sum.data = fread(paste0(temp.dir,"m_summary_chr_",j,"_rho_",l,".out.P",i_rep,".assoc.linear"))
  sum.data.infor = sum.data[,1:6]
  sum.data.list = list()
  for(i_rep in 1:(end-start+1)){
    sum.data = fread(paste0(temp.dir,"m_summary_chr_",j,"_rho_",l,".out.P",i_rep,".assoc.linear"))
    effect = sum.data[,c(7:9)]
    sum.data.list[[i_rep]] = effect
  }
  effect = bind_cols(sum.data.list)
  sum.data = cbind(sum.data.infor,effect)
  write.table(sum.data,file = paste0(cur.dir,"m_summary_chr_",j,"_rho_",l,"_sub_",i1),row.names = F,col.names = T,quote=F)
  
  
  

# for(i in 1:5){
#   fam <- data.frame(fread(paste0(cur.dir,eth[i],"/all_chr.tag.fam")))
#     for(j in 1:22){
#     write.table(fam, file = paste0(cur.dir,eth[i],"/chr",j,".tag.fam"),row.names = F,col.names = F,quote=F)
#   }
# }
# 
# 
