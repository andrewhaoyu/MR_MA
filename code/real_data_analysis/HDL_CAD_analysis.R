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
  select(CHR,SNP,BP,A1,TEST,NMISS,BETA,STAT,P)
write.table(assoc,file = "/data/zhangh24/MR_MA/result/real_data_analysis/HDL/HDL_assoc",quote=F,row.names = F,col.names=T)


#load the clumped SNPs
clump_SNP = as.data.frame(fread("/data/zhangh24/MR_MA/result/real_data_analysis/HDL/HDL_clump.clumped",header=T))
clump_SNP = clump_SNP[,3,drop=F]
colnames(clump_SNP) = "SNP"
clump_SNP_infor = left_join(clump_SNP,
                            HDL.KG,by="SNP")


clump_SNP_infor = clump_SNP_infor %>% 
  rename(beta_ex = beta,
         se_ex = se,
         P_ex = P) %>% 
  select(SNP,reference_allele,other_allele,CHR,BP,beta_ex,se_ex,P_ex,chr.pos)

#load CAD summary level statistics
CAD <- as.data.frame(fread("./data/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt",header=T))
colnames(CAD)[10] = "pvalue"
CAD = CAD %>% 
  rename(SNP=snptestid) %>% 
  mutate(chr.pos= paste0(chr,":",bp_hg19))
out_clump_SNP = inner_join(clump_SNP_infor,
                           CAD,by="chr.pos") %>% 
  rename(SNP=SNP.x) %>% 
  select(SNP,CHR,BP,reference_allele,other_allele,beta_ex,se_ex,P_ex,noneffect_allele,effect_allele,effect_allele_freq,logOR,se_gc,pvalue)

#align the SNP to the reference SNP
idx <- which(out_clump_SNP$reference_allele!=out_clump_SNP$effect_allele)
out_clump_SNP$logOR[idx] = -out_clump_SNP$logOR[idx]
allele_temp = out_clump_SNP$effect_allele[idx]
out_clump_SNP$effect_allele[idx] = out_clump_SNP$noneffect_allele[idx]
out_clump_SNP$noneffect_allele[idx] =allele_temp 

idx.order <- order(as.numeric(out_clump_SNP$CHR),
                   as.numeric(out_clump_SNP$BP))
out_clump_SNP = out_clump_SNP[idx.order,]

out_clump_SNP = out_clump_SNP[idx.order,] %>% 
  mutate(Gamma = out_clump_SNP$logOR,
         var_Gamma = out_clump_SNP$se_gc^2,
         gamma = out_clump_SNP$beta_ex,
         var_gamma = out_clump_SNP$se_ex^2)




write.table(out_clump_SNP,file = "/data/zhangh24/MR_MA/result/real_data_analysis/LDL/HDL_CAD.aligned",row.names = F,col.names = T,quote=F)
