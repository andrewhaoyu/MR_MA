#LDL analysis 
#first run LD clumping to get genome-wide significant SNPs
#second run the IVW and beta estimate model
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
         SNP = rsid) 


assoc = LDL %>% 
  select(CHR,SNP,BP,A1,TEST,NMISS,BETA,STAT,P)
write.table(assoc,file = "/data/zhangh24/MR_MA/result/realdata/LDL/LDL_assoc",quote=F,row.names = F,col.names=T)

pthr = 5E-08
r2thr = 0.1
kbpthr = 500
LD.clump.code <- paste0("/data/zhangh24/plink --bfile /gpfs/gsfs11/users/zhangh24/KG.plink/KG.all.chr --clump /data/zhangh24/MR_MA/result/realdata/LDL/LDL_assoc --clump-p1 ",pthr," --clump-r2 ",r2thr,"  --clump-kb ",kbpthr," --out /data/zhangh24/MR_MA/result/real_data_analysis/LDL/LDL_clump")
#run the code in terminal


#load the clumped SNPs
clump_SNP = as.data.frame(fread("/gpfs/gsfs11/users/zhangh24/MR_MA/result/real_data_analysis/LDL/LDL_clump.clumped",header=T))
clump_SNP = clump_SNP[,3,drop=F]
colnames(clump_SNP) = "rsid"
clump_SNP_infor = left_join(clump_SNP,
                            LDL,by="rsid")
clump_SNP_infor = clump_SNP_infor %>% 
  mutate(A2=toupper(A2)) %>% 
  rename(beta_ex = beta,
         se_ex = se,
         P_ex = P) %>% 
  select(SNP,A1,A2,CHR,BP,beta_ex,se_ex,P_ex)

#load CAD summary level statistics
CAD <- as.data.frame(fread("./data/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt",header=T))
colnames(CAD)[10] = "pvalue"
CAD = CAD %>% 
  rename(SNP=snptestid)
out_clump_SNP = inner_join(clump_SNP_infor,
                          CAD,by="SNP") %>% select(SNP,CHR,BP,A1,A2,beta_ex,se_ex,P_ex,noneffect_allele,effect_allele,effect_allele_freq,logOR,se_gc,pvalue)

#align the SNP to the reference SNP
idx <- which(out_clump_SNP$A1!=out_clump_SNP$noneffect_allele)
out_clump_SNP$logOR[idx] = -out_clump_SNP$logOR[idx]
#take out the reference allele and other allele since it's the same as exposure A1,A2

out_clump_SNP = out_clump_SNP %>% 
  select(SNP,CHR,BP,A1,A2,beta_ex,se_ex,P_ex,effect_allele_freq,logOR,se_gc,pvalue) 
Gamma = out_clump_SNP$logOR
var_Gamma = out_clump_SNP$se_gc^2
gamma = out_clump_SNP$beta_ex
var_gamma = out_clump_SNP$se_ex^2
IVW_s(Gamma,var_Gamma,
      gamma,var_gamma)




RatioExact = function(Gamma,var_Gamma,gamma,var_gamma,n){
  n = length(Gamma)
  ratio_est = Gamma/gamma
  for(k in 1:)
  n.simu <- 1000000
  z_Gamma <- rnorm(n.simu,mean =0,sd =sqrt((n-1)*var_Gamma))
  z_gamma <- rnorm(n.simu,mean = sqrt(n)*gamma,sd = sqrt((n-1)*var_gamma))
  #z_Gamma <- rnorm(n.simu,mean =0,sd =sqrt(var_Gamma))
  #z_gamma <- rnorm(n.simu,mean = gamma,sd = sqrt(var_gamma))
  #z_gamma <- rnorm(n.simu,mean = sqrt(n)*alpha_G,sd = sqrt((n-1)*var_gamma))
  
  true_distribution <- z_Gamma/z_gamma
  q_result <- quantile(true_distribution,c(0.025,0.975))
  cover = ifelse(ratio_est>=q_result[1]&
                   ratio_est<=q_result[2],1,0)
  ci_low <- ratio_est-q_result[2]
  ci_high <- ratio_est-q_result[1]
  return(c(cover,ci_low,ci_high))
}



Meta = function(coef_vec,var_vec){
  meta_var = (sum(1/var_vec))^-1
  meta_coef = meta_var*sum(coef_vec/var_vec)
  return(c(meta_coef,meta_var))
}
IVW_s = function(Gamma,var_Gamma,gamma,var_gamma){
  p <- length(Gamma)
  raio_vec = rep(0,p)
  ratio_var_vec = rep(0,p)
  
    raio_vec = Gamma/gamma
    ratio_var_vec =  var_Gamma/gamma^2+var_gamma*Gamma^2/gamma^4
  
  
  Meta_result = Meta(raio_vec,ratio_var_vec)
  ratio_ivw =   Meta_result[1]
  ratio_ivw_var = Meta_result[2]
  coef_low = ratio_ivw-1.96*sqrt(ratio_ivw_var)
  coef_high = ratio_ivw+1.96*sqrt(ratio_ivw_var)
  
  return(c(ratio_ivw,ratio_ivw_var,
           coef_low,coef_high))
}

