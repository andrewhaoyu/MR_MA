#LDL analysis 
#first run LD clumping to get genome-wide significant SNPs
#second run the IVW and beta estimate model
#update date: 012720
#down load data (http://mccarthy.well.ox.ac.uk/publications/2015/ENGAGE_1KG/LDL_Meta_ENGAGE_1000G.txt.gz)
#load LDL summary level statistics
#A1 is the effect allele and A2 is the noneffect allele in LDL analysis
setwd("/data/zhangh24/MR_MA")
library(data.table)
library(dplyr)
library(tidyr)
KG.SNP <- as.data.frame(fread("/data/zhangh24/KG.plink/KG.all.chr.bim",header=F))
colnames(KG.SNP) <- c("CHR","SNP","Nothing","BP","Allele1","Allele2")
KG.SNP = KG.SNP %>% 
  mutate(chr.pos = paste0(CHR,":",BP)) %>% select(SNP,chr.pos)


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
#clumping only takes 3 minutes
write.table(LD.clump.code,file = "/data/zhangh24/MR_MA/code/real_data_analysis/LD.clump_LDL.sh",quote = F,row.names = F,col.names = F)





#load the clumped SNPs
clump_SNP = as.data.frame(fread("/gpfs/gsfs11/users/zhangh24/MR_MA/result/real_data_analysis/LDL/LDL_clump.clumped",header=T))
clump_SNP = clump_SNP[,3,drop=F]
colnames(clump_SNP) = "SNP"
clump_SNP_infor = left_join(clump_SNP,
                            LDL.KG,by="SNP")


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
                          CAD,by="chr.pos") 

out_clump_SNP = out_clump_SNP%>% 
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
out_clump_SNP = out_clump_SNP[idx.order,] %>% 
  mutate(Gamma = out_clump_SNP$logOR,
         var_Gamma = out_clump_SNP$se_gc^2,
         gamma = out_clump_SNP$beta_ex,
         var_gamma = out_clump_SNP$se_ex^2)




write.table(out_clump_SNP,file = "/data/zhangh24/MR_MA/result/real_data_analysis/LDL/LDL_CAD.aligned",row.names = F,col.names = T,quote=F)


pcut <- c(5E-08)
l <- length(pcut)
IVW_c_result <- rep("c",l)
MRweight_result <- rep("c",l)
MR_egger_result <- rep("c",l)
MR_median_result <- rep("c",l)
MR_raps_result <- rep("c",l)
MR_presso_result <- rep("c",l)
n.snp <- rep(0,l)
keep.snp <- rep(0,l)
for(k in 1:l){
  pdx <- which(out_clump_SNP$P_ex<=pcut[k])
  out_clump_SNP_temp = out_clump_SNP[pdx,]
  Gamma = out_clump_SNP_temp$logOR
  var_Gamma = out_clump_SNP_temp$se_gc^2
  gamma = out_clump_SNP_temp$beta_ex
  var_gamma = out_clump_SNP_temp$se_ex^2
  n.snp[k] <- length(Gamma)
  num = 3
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                       gamma,var_gamma)
  IVW_c_result[k] <- paste0(round(IVW_c_temp[[1]],num)," (",round(IVW_c_temp[[4]],num),")")
  MR_weight_temp = MRWeight(Gamma,var_Gamma,
                            gamma,var_gamma)
  MRweight_result[k] <- paste0(round(MR_weight_temp[[1]],num)," (",round(MR_weight_temp[[4]],num),")")
  MRInputObject <- mr_input(bx = gamma,
                            bxse = sqrt(var_gamma),
                            by = Gamma,
                            byse = sqrt(var_Gamma))
  median_result <- mr_median(MRInputObject,
                             weighting = "weighted",
                             distribution = "normal",
                             alpha = 0.05,
                             iterations = 10000,
                             seed = 314159265)
  MR_median_result[k] <- paste0(round(median_result$Estimate,num)," (",round(median_result$StdError,num),")")
  egger_result <- mr_egger(MRInputObject,
                           robust = FALSE,
                           penalized = FALSE,
                           correl = FALSE,
                           distribution = "normal",
                           alpha = 0.05)
  MR_egger_result[k] <- paste0(round(egger_result$Estimate,num)," (",round(egger_result$StdError.Est,num),")")
  # raps_result <- mr.raps(data = data.frame(beta.exposure = gamma,
  #                                          beta.outcome = Gamma,
  #                                          se.exposure = sqrt(var_gamma),
  #                                          se.outcome = sqrt(var_Gamma)),
  #                                           )
  #MR_raps_result[k] <- paste0(round(raps_result$beta.hat,num)," (",round(raps_result$beta.se,num),")")
  summary.data = data.frame(E1_effect = gamma,
                            E1_se = sqrt(var_gamma),
                            
                            Y_effect = Gamma,
                            Y_se = sqrt(var_Gamma))
  presso_result <- mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = summary.data, NbDistribution = 1000,  SignifThreshold = 0.05)
  mr_presso_result = paste0(round(presso_result$`Main MR results`[1,3],num)," (",round(presso_result$`Main MR results`[1,4],num),")")
}


LDL.result.summary <- data.frame(IVW_c_result,MRweight_result,MR_egger_result,MR_median_result
                                 ,MR_raps_result,mr_presso_result)
write.csv(LDL.result.summary,file = "/data/MR_MA/result/real_data_analysis/LDL/LDL_clump.clumped")





#IVW estimate using summary level statistics
IVW_c = function(Gamma,var_Gamma,gamma,var_gamma){
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

MRLR <- function(Gamma,var_Gamma,gamma,var_gamma){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  
  #first step
  model1 = lm(Gamma~gamma-1)
  coef_est = coefficients(model1)
  W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  
  coef_best = sum(Gamma*gamma*W_vec)/sum(gamma^2*W_vec)
  sigma_est  = sum((Gamma-coef_est*gamma)^2)/(K-1)
  
  W_vec = 1/(var_Gamma+coef_est^2*var_gamma)
  xwx_iv = 1/sum(gamma^2*W_vec)
  
  var_coef_est = sigma_est*xwx_iv*t(gamma)%*%diag(W_vec)%*%diag(W_vec)%*%gamma*xwx_iv
  
  coef_low <- coef_est+qt(0.025,(K-1))*sqrt(var_coef_est)
  coef_high <- coef_est+qt(0.975,(K-1))*sqrt(var_coef_est)
  
  
  return(list(coef_est,coef_low,coef_high))
}









type <- rep("c",length(Gamma))
type[AR_result[[4]]] <- "include"
type[AR_result[[5]]] <- "remove"
library(ggplot2)
data <- data.frame(gamma,Gamma,type)
ggplot(data,aes(gamma,Gamma))+geom_point(aes(gamma,Gamma,color=type))


prop <- (var_gamma*Gamma^2/gamma^4)/(var_Gamma/gamma^2+var_gamma*Gamma^2/gamma^4)
boxplot(prop)
keep.ind <- AR_result[[4]]
IVW_s(Gamma[keep.ind],var_Gamma[keep.ind],
      gamma[keep.ind],var_gamma[keep.ind])
ARMethod(Gamma[keep.ind],
                      var_Gamma[keep.ind],
                      gamma[keep.ind],
                      var_gamma[keep.ind])
Gamma = Gamma[keep.ind]
var_Gamma = var_Gamma[keep.ind]
gamma = gamma[keep.ind]
var_gamma = var_gamma[keep.ind]


Gamma_all = Gamma
var_Gamma_all = var_Gamma
gamma_all = gamma
var_gamma_all = var_gamma



for(k in 1:79){
  print(k)
  Gamma = Gamma_all[c(1:k)]
  var_Gamma = var_Gamma_all[c(1:k)]
  gamma = gamma_all[c(1:k)]
  var_gamma = var_gamma_all[c(1:k)]
  print(ARMethod(Gamma,var_Gamma,
           gamma,var_gamma))
  
}


QuacForm <- function(Gamma,var_Gamma,gamma,var_gamma,beta_plug){
  K <- length(Gamma)
  #return the plug in test results
  return(sum((Gamma-beta_plug*gamma)^2/(var_Gamma+beta_plug^2*var_gamma))-qchisq(0.95,K))
}

ARMethod <- function(Gamma,var_Gamma,gamma,var_gamma){
  K <- length(Gamma)
  keep.ind <- c(1:K)
  
  beta_seq <- seq(-5,5,by=0.001)
  #initial run
  quan_result <- rep(0,length(beta_seq))
  for(k in 1:length(beta_seq)){
    quan_result[k] <- QuacForm(Gamma,var_Gamma,gamma,var_gamma,beta_seq[k])
  }
  #get the best estimate
  coef_best <- beta_seq[which.min(quan_result)]
  #test whether there are any special value
  beta_plug = coef_best
  
  Gamma_update = Gamma
  var_Gamma_update = var_Gamma
  gamma_update = gamma
  var_gamma_update = var_gamma
  p_test_value = pchisq((Gamma_update-beta_plug*gamma_update)^2/(var_Gamma_update+beta_plug^2*var_gamma_update),1,lower.tail = F)
  p_adjust_test_value = p.adjust(p_test_value,method="none")  
  
  keep_update <- keep.ind
  while(min(p_adjust_test_value)<=0.05){
    print(min(p_adjust_test_value))
    idx <- which.min(p_adjust_test_value)
    keep.ind <- keep.ind[-idx]
    Gamma_update = Gamma_update[-idx]
    var_Gamma_update = var_Gamma_update[-idx]
    gamma_update = gamma_update[-idx]
    var_gamma_update = var_gamma_update[-idx]
    quan_result <- rep(0,length(beta_seq))
    for(k in 1:length(beta_seq)){
      quan_result[k] <- QuacForm(Gamma_update,var_Gamma_update,gamma_update,var_gamma_update,beta_seq[k])
    }
    coef_best <- beta_seq[which.min(quan_result)]
    p_test_value = pchisq((Gamma_update-coef_best*gamma_update)^2/(var_Gamma_update+coef_best^2*var_gamma_update),1,lower.tail = F)
    p_adjust_test_value = p.adjust(p_test_value,method="none") 
  }
  coef_est = coef_best
  #get the confidence interval
 
  idx <- which(quan_result<=0)
  length(idx)
  beta_ci_range <- beta_seq[idx]
  coef_low <- min(beta_ci_range)
  coef_high <- max(beta_ci_range)
  remove.id <- c(1:K)[c(1:K)%in%keep.ind==F]
  
  Gamma_keep = Gamma[keep.ind]
  var_Gamma_keep = var_Gamma[keep.ind]
  gamma_keep = gamma[keep.ind]
  var_gamma_keep = var_gamma[keep.ind]
  
  var_coef_est = solve(t(gamma_keep)%*%solve(diag(var_Gamma_keep+coef_est^2*var_gamma_keep))%*%(gamma_keep))
  
  coef_high_update <- coef_est+1.96*sqrt(var_coef_est)
  coef_low_update <- coef_est-1.96*sqrt(var_coef_est)
  
  return(list(coef_est,coef_low,coef_high,
              keep.ind,remove.id,coef_low_update,coef_high_update
              ))
}

Meta = function(coef_vec,var_vec){
  meta_var = (sum(1/var_vec))^-1
  meta_coef = meta_var*sum(coef_vec/var_vec)
  return(c(meta_coef,meta_var))
}


#IVW estimate using summary level statistics
IVW_s = function(Gamma,var_Gamma,gamma,var_gamma){
  p <- length(Gamma)
  raio_vec = rep(0,p)
  ratio_var_vec = rep(0,p)
  
  raio_vec = Gamma/gamma
  ratio_var_vec =  var_Gamma/gamma^2
  
  
  Meta_result = Meta(raio_vec,ratio_var_vec)
  ratio_ivw =   Meta_result[1]
  ratio_ivw_var = Meta_result[2]
  coef_low = ratio_ivw-1.96*sqrt(ratio_ivw_var)
  coef_high = ratio_ivw+1.96*sqrt(ratio_ivw_var)
  return(c(ratio_ivw,ratio_ivw_var,
           coef_low,coef_high))
}

#IVW estimate using summary level statistics
IVW_c = function(Gamma,var_Gamma,gamma,var_gamma){
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



Gamma = as.numeric(est[[1]])
var_Gamma = as.numeric(est[[2]])
gamma = as.numeric(est[[3]])
var_gamma = as.numeric(est[[4]])











sig_SNPs <- data.frame(out_clump_SNP$SNP,
                       as.numeric(out_clump_SNP$CHR),
                       as.numeric(out_clump_SNP$BP),
                       as.numeric(out_clump_SNP$P),
                       stringsAsFactors = F)
colnames(sig_SNPs) <- c("SNP","CHR","position",
                        "p.value")
PositionPruning <- function(sig_SNPs){
  sig_SNPs_temp =sig_SNPs
  filter_result = NULL
  temp.ind = 1
  while(nrow(sig_SNPs_temp)!=0){
    
    idx = which.min(sig_SNPs_temp$p.value)
    filter_result = rbind(filter_result,sig_SNPs_temp[idx,])
    
    position.range <- 500*10^3
    filter_result_position = sig_SNPs_temp$position[idx]
    filter_CHR = sig_SNPs_temp$CHR[idx]
    idx.cut <- which((sig_SNPs_temp$position>=filter_result_position-position.range)&(sig_SNPs_temp$position<=filter_result_position+position.range)&(sig_SNPs_temp$CHR==filter_CHR ))
    
    
    sig_SNPs_temp = sig_SNPs_temp[-idx.cut,]
    temp.ind = temp.ind+1
  }
  return(filter_result)
}

out_clump_SNP_clean <- PositionPruning(sig_SNPs)[,1,drop=F]
out_clump_SNP_clean <- left_join(out_clump_SNP_clean,
                                 out_clump_SNP,
                                 by="SNP")

#order the SNPs by chr and position
idx.order <- order(out_clump_SNP_clean$CHR,
                   out_clump_SNP_clean$BP)
out_clump_SNP_clean = out_clump_SNP_clean[idx.order,]
