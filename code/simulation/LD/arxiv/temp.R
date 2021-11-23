library(tidyverse)
LD.snp = LD.snp %>% 
  separate(sep = ":",col = "SNP",into =  c("rsid","pos","Allele1","Allele2"),remove=F)

cau_snp = cau_snp_list[[1]]

jdx = which(abs(cau_snp$physical.pos-40715303)<=500000)
cau_snp[jdx,]
