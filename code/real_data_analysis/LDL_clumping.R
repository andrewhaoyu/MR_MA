#LDL analysis 
#first run LD clumping to get genome-wide significant SNPs
#second run the IVW and beta estimate model
#update date: 012720
#down load data (http://mccarthy.well.ox.ac.uk/publications/2015/ENGAGE_1KG/LDL_Meta_ENGAGE_1000G.txt.gz)
#load LDL summary level statistics
#A1 is the effect allele and A2 is the noneffect allele in LDL analysis
#read KG SNP information
KG.SNP <- as.data.frame(fread("/data/zhangh24/KG.plink/EUR/chr_all.bim",header=F))
colnames(KG.SNP) <- c("CHR","SNP","Nothing","BP","Allele1","Allele2")
KG.SNP = KG.SNP %>% 
  mutate(chr.pos = paste0(CHR,":",BP)) %>% select(SNP,chr.pos)





LDL = as.data.frame(fread("./data/LDL_Meta_ENGAGE_1000G.txt",header=T))
n = nrow(LDL)
colnames(LDL)[7] = "P"
LDL = LDL %>% 
  mutate(chr = gsub("chr","",chr)) %>%mutate(chr.pos = paste0(chr,":",pos))%>% 
  mutate(A1=toupper(reference_allele),
         TEST = "ADD",
         NMISS= 0,
         BETA=beta,
         STAT = rnorm(n)) %>% rename(CHR=chr,
                                     BP =pos)
#get the SNPs that are shared by LDL and KG
LDL.KG = inner_join(LDL,KG.SNP,
                    by="chr.pos")

assoc = LDL.KG %>% 
  select(CHR,SNP,BP,A1,TEST,NMISS,BETA,STAT,P)
write.table(assoc,file = "/data/zhangh24/MR_MA/result/real_data_analysis/LDL/LDL_assoc",quote=F,row.names = F,col.names=T)

pthr = 5E-08
r2thr = 0.001
kbpthr = 3000
LD.clump.code <- paste0("/data/zhangh24/software/plink2 --bfile /data/zhangh24/KG.plink/EUR/chr_all --clump /data/zhangh24/MR_MA/result/real_data_analysis/LDL/LDL_assoc --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /data/zhangh24/MR_MA/result/real_data_analysis/LDL/LDL_clump")
#run the code in terminal
write.table(LD.clump.code,file = "/data/zhangh24/MR_MA/code/real_data_analysis/LD.clump_LDL.sh",quote = F,row.names = F,col.names = F)
