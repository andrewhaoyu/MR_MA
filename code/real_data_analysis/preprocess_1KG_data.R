#Goal: preprocess KG information
#read KG SNP information
KG.SNP <- as.data.frame(fread("/data/zhangh24/KG.plink/EUR/chr_all.bim",header=F))

colnames(KG.SNP) <- c("CHR","SNP","Nothing","BP","Allele1","Allele2")
#load KG SNPs MAF information
library(data.table)
library(dplyr)
snp.infor.list = list()
for(i in 1:22){
  leg <- as.data.frame(fread(paste0("/data/zhangh24/KG.impute2/1000GP_Phase3/1000GP_Phase3_chr",i,".legend"),header=T))
  snp.infor.list[[i]] = leg %>% 
    rename(SNP=id,EAF=EUR) %>% 
    select(SNP,EAF)
  
  
}
snp.infor.select = rbindlist(snp.infor.list)



KG.SNP = left_join(KG.SNP,snp.infor.select,by="SNP")

KG.SNP = KG.SNP %>% 
  mutate(chr.pos = paste0(CHR,":",BP)) %>% select(SNP,chr.pos,EAF) 
#save KG information
save(KG.SNP, file = "/data/zhangh24/KG.plink/EUR/KG.SNP_process.Rdata")
