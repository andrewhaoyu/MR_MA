setwd("/data/zhangh24/MR_MA")
library(data.table)
library(dplyr)
library(tidyr)

#read KG SNP information
KG.SNP <- as.data.frame(fread("/data/zhangh24/KG.plink/EUR/chr_all.bim",header=F))
colnames(KG.SNP) <- c("CHR","SNP","Nothing","BP","Allele1","Allele2")
KG.SNP = KG.SNP %>% 
  mutate(chr.pos = paste0(CHR,":",BP)) %>% select(SNP,chr.pos)


bmi = as.data.frame(fread("./data/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt",header=T))
n = nrow(bmi)

bmi = bmi %>% 
  mutate(chr.pos = paste0(CHR,":",POS))%>% 
  mutate(TEST = "ADD",
         NMISS= 0,
         STAT = BETA/SE)
#get the SNPs that are shared by bmi and KG
bmi.KG = inner_join(bmi,KG.SNP,
                    by="chr.pos")



bmi.KG =  bmi.KG %>% 
  rename(SNP=SNP.y,
         BP = POS,
         A1 = Tested_Allele)


assoc = bmi.KG %>% 
  select(CHR,SNP,BP,A1,TEST,NMISS,BETA,STAT,P) %>% 
  filter(P<=5E-08)
#idx <- which(assoc$P<=5E-08)
write.table(assoc,file = "/data/zhangh24/MR_MA/result/real_data_analysis/bmi/bmi_assoc",quote=F,row.names = F,col.names=T)

pthr = 0.1
r2thr = 0.001
kbpthr = 3000

  
LD.clump.code <- paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/EUR/chr_all --clump /data/zhangh24/MR_MA/result/real_data_analysis/bmi/bmi_assoc --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /data/zhangh24/MR_MA/result/real_data_analysis/bmi/bmi_clump")
#run the code in terminal
write.table(LD.clump.code,file = "/data/zhangh24/MR_MA/code/real_data_analysis/LD.clump_bmi.sh",quote = F,row.names = F,col.names = F)
