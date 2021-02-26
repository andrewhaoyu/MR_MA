#HDL analysis 
#first run LD clumping to get genome-wide significant SNPs
#second run the IVW and beta estimate model
#update date: 012720
#down load data (http://mccarthy.well.ox.ac.uk/publications/2015/ENGAGE_1KG/HDL_Meta_ENGAGE_1000G.txt.gz)
#load HDL summary level statistics
#A1 is the effect allele and A2 is the noneffect allele in HDL analysis
setwd("/data/zhangh24/MR_MA")
library(data.table)
library(dplyr)
library(tidyr)

#read KG SNP information
KG.SNP <- as.data.frame(fread("/data/zhangh24/KG.plink/EUR/chr_all.bim",header=F))
colnames(KG.SNP) <- c("CHR","SNP","Nothing","BP","Allele1","Allele2")
KG.SNP = KG.SNP %>% 
  mutate(chr.pos = paste0(CHR,":",BP)) %>% select(SNP,chr.pos)



HDL = as.data.frame(fread("./data/HDL_Meta_ENGAGE_1000G.txt",header=T))
n = nrow(HDL)
colnames(HDL)[7] = "P"
HDL = HDL %>% 
  mutate(chr = gsub("chr","",chr)) %>%mutate(chr.pos = paste0(chr,":",pos))%>% 
  mutate(A1=toupper(reference_allele),
         TEST = "ADD",
         NMISS= 0,
         BETA=beta,
         STAT = rnorm(n)) %>% rename(CHR=chr,
                                     BP =pos)
#get the SNPs that are shared by HDL and KG
HDL.KG = inner_join(HDL,KG.SNP,
                    by="chr.pos")


assoc = HDL.KG %>% 
  select(CHR,SNP,BP,A1,TEST,NMISS,BETA,STAT,P) %>% 
  filter(P<=5E-08)
write.table(assoc,file = "/data/zhangh24/MR_MA/result/real_data_analysis/HDL/HDL_assoc",quote=F,row.names = F,col.names=T)

pthr = 5E-08
r2thr = 0.001
kbpthr = 3000
LD.clump.code <- paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/EUR/chr_all --clump /data/zhangh24/MR_MA/result/real_data_analysis/HDL/HDL_assoc --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /data/zhangh24/MR_MA/result/real_data_analysis/HDL/HDL_clump")
#run the code in terminal
write.table(LD.clump.code,file = "/data/zhangh24/MR_MA/code/real_data_analysis/LD.clump_HDL.sh",quote = F,row.names = F,col.names = F)
