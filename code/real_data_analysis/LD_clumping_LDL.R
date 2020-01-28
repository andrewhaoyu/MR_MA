#Merge the ratio estiamtes results
#update date: 012720
#down load data (http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_LDL.txt.gz)
#load LDL summary level statistics
setwd("/gpfs/gsfs11/users/zhangh24/MR_MA")
library(data.table)
LDL = as.data.frame(fread("./data/jointGwasMc_LDL.txt",header=T))
library(dplyr)
library(tidyr)
n = nrow(LDL)
colnames(LDL)[9] = "P"
LDL = LDL %>% 
  mutate(chr.pos = gsub("chr","",LDL$SNP_hg19)) %>% 
  separate(chr.pos,c("CHR","BP"),":") %>% 
  mutate(A1=toupper(A1),
         TEST = "ADD",
         NMISS= 0,
         BETA=beta,
         STAT = rnorm(n),
         SNP = rsid) %>% 
  select(CHR,SNP,BP,A1,TEST,NMISS,BETA,STAT,P)
write.table(LDL,file = "/data/zhangh24/MR_MA/result/realdata/LDL/LDL_assoc",quote=F,row.names = F,col.names=T)

pthr = 5E-08
r2thr = 0.1
kbpthr = 500
LD.clump.code <- paste0("/data/zhangh24/plink --bfile /gpfs/gsfs11/users/zhangh24/KG.plink/ --clump /data/zhangh24/MR_MA/result/realdata/LDL/LDL_assoc --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /data/zhangh24/MR_MA/result/realdata/LDL/LDL_clump")
