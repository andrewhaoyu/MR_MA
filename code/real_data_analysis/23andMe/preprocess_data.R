#Function for transforming fasrc array to vec
TransArraytoVec <- function(array,vec){
  n = length(vec)
  command = paste0("expand.grid(")
  for(k in 1:(n-1)){
    temp = paste0("c(1:",vec[k],"),")
    command = paste0(command,temp)
  }
  temp = paste0("c(1:",vec[n],"))")
  command = paste0(command,temp)
 result =  eval(parse(text = command))
  return(result[array,])
}
setwd("/n/holystore01/LABS/xlin/Lab/hzhang/MR_MA/")
library(data.table)
library(R.utils)
library(tidyverse)
args = commandArgs(trailingOnly = T)
array = as.numeric(args[[1]])
vecargs = TransArraytoVec(array,c(4,5))
i = as.numeric(vecargs[1])
l = as.numeric(vecargs[2])
eth = c("european","african_american",
        "latino","trans_ethnic")
traits = c("positive_vs_negative",
           "positive_hospitalized_dx_negative_controls",
           "positive_pneumonia_dx_negative_controls",
           "positive_respiratory_broad_dx_negative_controls",
           "positive_respiratory_support_dx_negative_controls")
snp.infor <- fread("./data/covid/8.2_Annotation/all_snp_info.txt")
data = fread(paste0("./data/covid/",eth[i],"/stats/covid_test_",traits[l],".dat.gz"))
imputed.infor = fread("./data/covid/8.2_Annotation/im_snp_stat.txt")
data.pass = data %>% filter(pass=="Y") %>% 
  select(all.data.id,src,pvalue,effect,stderr)
snp.infor = snp.infor %>% 
  select(all.data.id,im.data.id,assay.name,scaffold,position,alleles)
imputed.infor = imputed.infor %>% 
  select(im.data.id,freq.a,freq.b)
snp.infor.select = left_join(snp.infor,imputed.infor,by="im.data.id")
data.com = left_join(data.pass,snp.infor.select,by="all.data.id")

data.com = data.com %>% 
  mutate(CHR = gsub("chr","",scaffold),
         MAF = ifelse(freq.b<=0.5,freq.b,1-freq.b)) %>% 
  filter(MAF>=0.01) %>% 
  mutate(CHR = ifelse(CHR=="X",23,CHR))

data.com = data.com %>% 
  rename(rsid = assay.name) %>% 
  separate(alleles,into = c("non_effect_allele","effect_allele"),sep = "/")
 
data = data.com %>% 
  select(all.data.id,pvalue,effect,stderr,im.data.id,
         rsid,position,non_effect_allele,effect_allele,freq.a,freq.b,CHR,MAF) %>% 
  rename(effect_allele_freq = freq.b,
         non_effect_allele_freq = freq.a)

save(data,file = paste0("./data/cleaned/",eth[i],"/",traits[l],".rdata"))
