#goal subset the SNPs data to hm3 list
i =1
j = 22
library(data.table)
library(dplyr)

sid<-Sys.getenv('SLURM_JOB_ID')
dir.create(paste0('/lscratch/',sid,'/test'),showWarnings = FALSE)
temp.dir = paste0('/lscratch/',sid,'/test/')
eth <- c("EUR","AFR","AMR","EAS","SAS")
cur.dir <- "/data/zhangh24/multi_ethnic/result/LD_simulation_new/"
system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.bed ",temp.dir,eth[i],"chr",j,".tag.bed"))
system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.bim ",temp.dir,eth[i],"chr",j,".tag.bim"))
system(paste0("cp ",cur.dir,eth[i],"/chr",j,".tag.fam ",temp.dir,eth[i],"chr",j,".tag.fam"))

MAF.cutoff = 0.005
if(i==1){
  MAF.cutoff = 0.01
}
library(rlang)
setwd("/data/zhangh24/multi_ethnic/")
load("./result/LD_simulation_new/snp.infor.match37_38.rdata")
snp.infor = snp.infor.match %>% 
  filter(!!sym(eth[i])>=MAF.cutoff&
           !!sym(eth[i])<=(1-MAF.cutoff)&
           CHR==j) %>% 
  rename(SNP=id) %>% 
  select(SNP,rs_id,EUR)

hm3.list <- as.data.frame(fread(paste0(cur.dir,"hm3rsid.txt"),header  =T))
#mega.infor <- as.data.frame(fread(paste0(cur.dir,"snpBatch_ILLUMINA_1062317")))
#colnames(mega.infor)[5] <- "rsid"
colnames(hm3.list)[1] = "rs_id"

snp.infor.subset = inner_join(snp.infor,hm3.list,by="rs_id") %>% 
  select(SNP,EUR)
#save(snp.infor.subset,file = "/data/zhangh24/MR_MA/result/LD/chr22_snp_infor.rdata")
write.table(snp.infor.subset,file = paste0(temp.dir,"extract_snp_list.txt"),row.names = F,col.names = F,quote=F)

snp.infor.subset = inner_join(snp.infor,hm3.list,by="rs_id") %>% 
  mutate(MAF = ifelse(EUR<0.5,EUR,1-EUR)) %>% 
  select(SNP,rs_id,MAF)
save(snp.infor.subset,file ="/data/zhangh24/MR_MA/result/LD/chr22_snp_infor.rdata")
all.fam <- as.data.frame(fread(paste0(temp.dir,eth[i],"chr",j,".tag.fam")))
sub.fam <- all.fam[1:120000,]
write.table(sub.fam,file = paste0(temp.dir,"sub_fam.txt"),row.names = F,col.names = F,quote=F)
res = system(paste0("/data/zhangh24/software/plink2 --bfile ",temp.dir,eth[i],"chr",j,".tag --extract ",temp.dir,"extract_snp_list.txt --out ",temp.dir,"chr",j,".hm3 --make-bed --keep ",temp.dir,"sub_fam.txt"))

if(res==2){
  stop()
}
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
system(paste0("mv ",temp.dir,"/chr",j,".hm3.bed ",cur.dir,"chr",j,".hm3.bed"))
system(paste0("mv ",temp.dir,"/chr",j,".hm3.bim ",cur.dir,"chr",j,".hm3.bim"))
system(paste0("mv ",temp.dir,"/chr",j,".hm3.fam ",cur.dir,"chr",j,".hm3.fam"))




#subset 3000 people for clumping purpose
setwd("/data/zhangh24/MR_MA/result/LD")
cur.dir = "/data/zhangh24/MR_MA/result/LD/"
j = 22
all.fam <- as.data.frame(fread(paste0("chr",j,".hm3.fam")))
idx <- sample(c(1:12000),3000)
sub.fam <- all.fam[idx,]
write.table(sub.fam,file = paste0("sub_fam.txt"),row.names = F,col.names = F,quote=F)
res = system(paste0("/data/zhangh24/software/plink2 --bfile ",cur.dir,"chr",j,".hm3  --out ",cur.dir,"chr",j,".sub.hm3 --make-bed --keep ",cur.dir,"sub_fam.txt"))
