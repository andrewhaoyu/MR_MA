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



#MR-analysis for BMI-BC
#load the clumped SNPs
clump_SNP = as.data.frame(fread("/data/zhangh24/MR_MA/result/real_data_analysis/bmi/bmi_clump.clumped",header=T))
clump_SNP = clump_SNP[,3,drop=F]
colnames(clump_SNP) = "SNP"
clump_SNP_infor = left_join(clump_SNP,
                            bmi.KG,by="SNP")


clump_SNP_infor = clump_SNP_infor %>% 
  rename(beta_ex = BETA,
         se_ex = SE,
         P_ex = P,
         reference_allele = A1,
         other_allele = Other_Allele) %>% 
  select(SNP,reference_allele,other_allele,CHR,BP,beta_ex,se_ex,P_ex,chr.pos)

#load BC summary level statistics
BC<- as.data.frame(fread("./data/icogs_onco_gwas_meta_overall_breast_cancer_summary_level_statistics.txt",header=T))

BC = BC %>% 
  mutate(chr.pos= paste0(chr.Onco,":",Position.Onco))
out_clump_SNP = inner_join(clump_SNP_infor,
                           BC,by="chr.pos") %>% 
  select(SNP,CHR,BP,reference_allele,other_allele,beta_ex,se_ex,P_ex,Baseline.Onco,Effect.Meta,Baseline.Meta,Beta.meta,sdE.meta,p.meta)

#align the SNP to the reference SNP
idx <- which(out_clump_SNP$reference_allele!=out_clump_SNP$Effect.Meta)
out_clump_SNP$Beta.meta[idx] = -out_clump_SNP$Beta.meta[idx]
allele_temp = out_clump_SNP$Effect.Meta[idx]
out_clump_SNP$Effect.Meta[idx] = out_clump_SNP$Baseline.Meta[idx]
out_clump_SNP$Baseline.Meta[idx] =allele_temp 

idx.order <- order(as.numeric(out_clump_SNP$CHR),
                   as.numeric(out_clump_SNP$BP))
out_clump_SNP = out_clump_SNP[idx.order,]
out_clump_SNP = out_clump_SNP[idx.order,] %>% 
  mutate(Gamma = out_clump_SNP$Beta.meta,
         var_Gamma = out_clump_SNP$sdE.meta^2,
         gamma = out_clump_SNP$beta_ex,
         var_gamma = out_clump_SNP$se_ex^2)




write.table(out_clump_SNP,file = "/data/zhangh24/MR_MA/result/real_data_analysis/LDL/BMI_BC.aligned",row.names = F,col.names = T,quote=F)




#load breast_cancer_subtypes_data
BC = as.data.frame(fread(paste0("/data/zhangh24/breast_cancer_data_analysis/discovery_SNP/prepare_summary_level_statistics/result/icogs_onco_meta_intrinsic_subtypes_summary_level_statistics.txt")))
BC = BC %>% 
  mutate(chr.pos= paste0(chr.Onco,":",Position.Onco))
out_clump_SNP = inner_join(clump_SNP_infor,
                           BC,by="chr.pos") %>% 
  select(SNP,CHR,BP,reference_allele,other_allele,beta_ex,se_ex,P_ex,
         Baseline.Onco,
         Effect.Meta,
         Baseline.Meta,
         luminal_A_log_or_meta,
         Luminal_A_sd_meta,
         Luminal_B_log_or_meta,
         Luminal_B_sd_meta,
         Luminal_B_HER2Neg_log_or_meta,
         Luminal_B_HER2Neg_sd_meta,
         HER2_Enriched_log_or_meta,
         HER2_Enriched_sd_meta,
         Triple_Neg_log_or_meta,
         Triple_Neg_sd_meta)


colSums(is.na(out_clump_SNP))

out_clump_SNP = out_clump_SNP %>%
  mutate(luminal_A_log_or_meta=ifelse(reference_allele==Effect.Meta,luminal_A_log_or_meta,-luminal_A_log_or_meta),
         Luminal_B_log_or_meta=ifelse(reference_allele==Effect.Meta,Luminal_B_log_or_meta,-Luminal_B_log_or_meta),
         Luminal_B_HER2Neg_log_or_meta =ifelse(reference_allele==Effect.Meta,Luminal_B_HER2Neg_log_or_meta,-Luminal_B_HER2Neg_log_or_meta),
         HER2_Enriched_log_or_meta = ifelse(reference_allele==Effect.Meta,HER2_Enriched_log_or_meta,-HER2_Enriched_log_or_meta),
         Triple_Neg_log_or_meta = ifelse(reference_allele==Effect.Meta,Triple_Neg_log_or_meta,-Triple_Neg_log_or_meta))

idx.order <- order(as.numeric(out_clump_SNP$CHR),
                   as.numeric(out_clump_SNP$BP))
out_clump_SNP = out_clump_SNP[idx.order,]




write.table(out_clump_SNP,file = "/data/zhangh24/MR_MA/result/real_data_analysis/bmi/BMI_BCsub.aligned",row.names = F,col.names = T,quote=F)





