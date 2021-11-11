args = commandArgs(trailingOnly = T)
i = as.numeric(args[[1]])
l = as.numeric(args[[2]])
#i1 represent the sub id, split 100 replciates into 10
i1 = as.numeric(args[[3]])
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
#use a subset of 3000 people for clumping purpose
system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.bed /lscratch/",sid,"/test/chr",j,".sub.hm3.bed"))
system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.bim /lscratch/",sid,"/test/chr",j,".sub.hm3.bim"))
system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.fam /lscratch/",sid,"/test/chr",j,".sub.hm3.fam"))
sum.data.m = as.data.frame(fread(paste0(cur.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l)))
pthr = 0.1
r2thr = 0.01
kbpthr = 500
for(i_rep in  start:end){
  p = sum.data.m[,(6+3*i_rep)]
  sum.data.sub = sum.data.m[,c(1:6,6+(i_rep-1)*3+c(1:3))]
  colnames(sum.data.sub)[7:9] =  c("BETA","STAT","P")
 
  write.table(sum.data.sub,file =paste0(temp.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_rep_",i_rep),col.names = T,row.names = F,quote=F)
 
  res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,"chr",j,".hm3 --clump ",temp.dir,"m_summary_chr_",j,"_rho_",l,"_rep_",i_rep," --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"LD_chr_",j,"_rho_",l,"_rep_",i_rep))
  clump.result = as.data.frame(fread(paste0(temp.dir,"LD_chr_",j,"_rho_",l,"_rep_",i_rep,".clumped")))
  clump.snp = clump.result[,3,drop=F]
  write.table(clump.snp,file = paste0(cur.dir,"LD_chr_",j,"_rho_",l,"_rep_",i_rep,".clumped"),row.names = F,col.names = T,quote=F)
}
  