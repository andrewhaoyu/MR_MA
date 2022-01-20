args = commandArgs(trailingOnly = T)
i = 1
l = as.numeric(args[[1]])
r_ind = as.numeric(args[[2]])
#i1 represent the sub id, split 100 replciates into 10
# i1 = as.numeric(args[[3]])
j = 22
num = 1
library(dplyr)
library(bc2)
library(data.table)
n.rep = n_rep = 100
#start.end = startend(n.rep,num,i1)
start = 1
end = 100
#load the phenotpypes data and use plink to run
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
#use a subset of 3000 people for clumping purpose
system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.bed /lscratch/",sid,"/test/chr",j,".sub.hm3.bed"))
system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.bim /lscratch/",sid,"/test/chr",j,".sub.hm3.bim"))
system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.fam /lscratch/",sid,"/test/chr",j,".sub.hm3.fam"))

r2_vec = c(0.001,0.2,0.4,0.6)
r2 = r2_vec[r_ind]
# for(i in 1:2){
#   for(l in 1:3){
#for(v in 1:3){
 v = 1 
  
  sum.data.m = as.data.frame(fread(paste0(cur.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v)))
  pthr = 0.1
  r2thr = r2
  kbpthr = 500
  for(i_rep in  start:end){
    p = sum.data.m[,(6+3*i_rep)]
    sum.data.sub = sum.data.m[,c(1,3,2,6,6+(i_rep-1)*3+c(1:3))]
    colnames(sum.data.sub) =  c("CHR","SNP","BP","A1","BETA","SE","P")
    write.table(sum.data.sub,file =paste0(temp.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep),col.names = T,row.names = F,quote=F)
    
    res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,"chr",j,".sub.hm3 --clump ",temp.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep," --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"LD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep))
    clump.result = as.data.frame(fread(paste0(temp.dir,"LD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep,".clumped")))
    clump.snp = clump.result[,3,drop=F]
    write.table(clump.snp,file = paste0(cur.dir,"LD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep,"_rvec_",r_ind,".clumped"),row.names = F,col.names = T,quote=F)
  }
#}
#   }
# }


# 
# j = 22
# num = 1
# library(dplyr)
# library(bc2)
# library(data.table)
# n.rep = n_rep = 100
# #start.end = startend(n.rep,num,i1)
# start = 1
# end = 100
# #load the phenotpypes data and use plink to run
# cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
# sid<-Sys.getenv('SLURM_JOB_ID')
# dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
# temp.dir = paste0('/lscratch/',sid,'/test/')
# #use a subset of 3000 people for clumping purpose
# system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.bed /lscratch/",sid,"/test/chr",j,".sub.hm3.bed"))
# system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.bim /lscratch/",sid,"/test/chr",j,".sub.hm3.bim"))
# system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.fam /lscratch/",sid,"/test/chr",j,".sub.hm3.fam"))
# # for(i in 1:2){
# #   for(l in 1:3){
# for(v in 1:3){
#   
#   
#   sum.data.m = as.data.frame(fread(paste0(cur.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v)))
#   pthr = 0.1
#   r2thr = 0.6
#   kbpthr = 500
#   for(i_rep in  start:end){
#     p = sum.data.m[,(6+3*i_rep)]
#     sum.data.sub = sum.data.m[,c(1,3,2,6,6+(i_rep-1)*3+c(1:3))]
#     colnames(sum.data.sub) =  c("CHR","SNP","BP","A1","BETA","SE","P")
#     write.table(sum.data.sub,file =paste0(temp.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep),col.names = T,row.names = F,quote=F)
#     
#     res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,"chr",j,".sub.hm3 --clump ",temp.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep," --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"LD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep))
#     clump.result = as.data.frame(fread(paste0(temp.dir,"LD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep,".clumped")))
#     clump.snp = clump.result[,3,drop=F]
#     write.table(clump.snp,file = paste0(cur.dir,"strongLD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep,".clumped"),row.names = F,col.names = T,quote=F)
#   }
# }
#   }
# }



# 
# j = 22
# num = 1
# library(dplyr)
# library(bc2)
# library(data.table)
# n.rep = n_rep = 100
# #start.end = startend(n.rep,num,i1)
# start = 1
# end = 100
# #load the phenotpypes data and use plink to run
# cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
# sid<-Sys.getenv('SLURM_JOB_ID')
# dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
# temp.dir = paste0('/lscratch/',sid,'/test/')
# #use a subset of 3000 people for clumping purpose
# system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.bed /lscratch/",sid,"/test/chr",j,".sub.hm3.bed"))
# system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.bim /lscratch/",sid,"/test/chr",j,".sub.hm3.bim"))
# system(paste0("cp ", cur.dir,"chr",j,".sub.hm3.fam /lscratch/",sid,"/test/chr",j,".sub.hm3.fam"))
# # for(i in 1:2){
# #   for(l in 1:3){
# for(v in 1:3){
#   
#   
#   sum.data.m = as.data.frame(fread(paste0(cur.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v)))
#   pthr = 0.1
#   r2thr = 0.3
#   kbpthr = 500
#   for(i_rep in  start:end){
#     p = sum.data.m[,(6+3*i_rep)]
#     sum.data.sub = sum.data.m[,c(1,3,2,6,6+(i_rep-1)*3+c(1:3))]
#     colnames(sum.data.sub) =  c("CHR","SNP","BP","A1","BETA","SE","P")
#     write.table(sum.data.sub,file =paste0(temp.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep),col.names = T,row.names = F,quote=F)
#     
#     res = system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile ",temp.dir,"chr",j,".sub.hm3 --clump ",temp.dir,"m_summary_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep," --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out ",temp.dir,"LD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep))
#     clump.result = as.data.frame(fread(paste0(temp.dir,"LD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep,".clumped")))
#     clump.snp = clump.result[,3,drop=F]
#     write.table(clump.snp,file = paste0(cur.dir,"medLD_chr_",j,"beta_",i,"_rho_",l,"_ple_",v,"_rep_",i_rep,".clumped"),row.names = F,col.names = T,quote=F)
#   }
# }
#   }
# }
