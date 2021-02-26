pcut <- c(5E-08)
trait_vec = c("LDL-CAD",
              "HDL-CAD",
              "BMI-CAD",
              "LDL-BC",
              "HDL-BC",
              "BMI-BC"
)
total <- length(trait_vec)
IVW_c_result <- rep("c",total)

MRweight_result <- rep("c",total)
MR_egger_result <- rep("c",total)
MR_median_result <- rep("c",total)
MR_raps_result <- rep("c",total)
MR_presso_result <- rep("c",total)
n.snp <- rep(0,total)
keep.snp <- rep(0,total)

files = c("./result/real_data_analysis/LDL_CAD.aligned",
             "./result/real_data_analysis/HDL_CAD.aligned",
          "./result/real_data_analysis/BMI_CAD.aligned",
          "./result/real_data_analysis/LDL_BC.aligned",
          "./result/real_data_analysis/HDL_BC.aligned",
          "./result/real_data_analysis/BMI_BC.aligned")
for(k in 4:6){
  filename = files[k]
  out_clump_SNP = read.table(filename,header=T)
  #out_clump_SNP_temp = out_clump_SNP[pdx,]
  Gamma = out_clump_SNP$Gamma
  var_Gamma = out_clump_SNP$var_Gamma
  gamma = out_clump_SNP$gamma
  var_gamma = out_clump_SNP$var_gamma
  
  num = 3
  IVW_c_temp <- IVW_c(Gamma,var_Gamma,
                      gamma,var_gamma)
  IVW_c_result[k] <- paste0(round(IVW_c_temp[[1]],num)," (",round(IVW_c_temp[[4]],num),")")
  MR_weight_temp = MRWeight(Gamma,var_Gamma,
                            gamma,var_gamma)
  keep.id = MR_weight_temp[[5]]
  out.id = MR_weight_temp[[6]]
  
  keep.ind = rep(FALSE,length(Gamma))
  keep.ind[out.id] = TRUE
  keep.ind = factor(keep.ind,levels=c(TRUE,FALSE))
  lm(Gamma[keep.id]~gamma[keep.id])
  
  beta_est = MR_weight_temp[[1]]
  plot.data = data.frame(gamma,Gamma,Weight = 1/(var_Gamma+beta_est^2*var_gamma),Removed = keep.ind)
  p = ggplot(plot.data,aes(gamma,Gamma))+
    geom_point(aes(col=Removed,size = Weight),alpha = 0.9)+
    geom_abline(slope = MR_weight_temp[[1]],intercept = 0)+
    theme_Publication()+
    scale_colour_Publication()+
    ggtitle(paste0(trait_vec[k]))+
    xlab("Alpha")+
    geom_abline(slope =-0.502,intercept = 0.006, linetype = "dashed")
  
  png(filename = paste0("./result/real_data_analysis/",trait_vec[k],"annotate.png"),width = 8,height = 6, units = "in",res= 300)
  print(p)
  dev.off()
  
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
  raps_result <- mr.raps(data = data.frame(beta.exposure = gamma,
                                           beta.outcome = Gamma,
                                           se.exposure = sqrt(var_gamma),
                                           se.outcome = sqrt(var_Gamma)),
                                            )
  MR_raps_result[k] <- paste0(round(raps_result$beta.hat,num)," (",round(raps_result$beta.se,num),")")
  summary.data = data.frame(E1_effect = gamma,
                            E1_se = sqrt(var_gamma),
                            
                            Y_effect = Gamma,
                            Y_se = sqrt(var_Gamma))
  presso_result <- mr_presso(BetaOutcome = "Y_effect", BetaExposure = "E1_effect", SdOutcome = "Y_se", SdExposure = "E1_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = summary.data, NbDistribution = 1000,  SignifThreshold = 0.05)
  mr_presso_result[k] = paste0(round(presso_result$`Main MR results`[1,3],num)," (",round(presso_result$`Main MR results`[1,4],num),")")
}
result = data.frame(IVW  = IVW_c_result,
                    MR_Egger = MR_egger_result,
                    MR_Median = MR_median_result,
                    MR_Presso = mr_presso_result,
                    RAPs = MR_raps_result,
                    MR_Weight = MRweight_result)
write.csv(result,file = "./result/real_data_analysis/all_real_data_analysis.csv")




#analysis only for MR-weight

pcut <- c(5E-08)
trait_vec = c("LDL-CAD",
              "HDL-CAD",
              "BMI-CAD",
              "LDL-BC",
              "HDL-BC",
              "BMI-BC",
              "LDL-Luminal_A",
              "LDL-Luminal_B",
              "LDL-Luminal_B_HER2Neg",
              "LDL-HER2_Enriched",
              "LDL-TN",
              "HDL-Luminal_A",
              "HDL-Luminal_B",
              "HDL-Luminal_B_HER2Neg",
              "HDL-HER2_Enriched",
              "HDL-TN",
              "BMI-Luminal_A",
              "BMI-Luminal_B",
              "BMI-Luminal_B_HER2Neg",
              "BMI-HER2_Enriched",
              "BMI-TN"
)

total <- length(trait_vec)

MRweight_result <- rep(0,total)
MRweight_p_value <- rep(0,total)
Exposure = rep("c",total)
Outcome = rep("c",total)

files = c("./result/real_data_analysis/LDL_CAD.aligned",
          "./result/real_data_analysis/HDL_CAD.aligned",
          "./result/real_data_analysis/BMI_CAD.aligned",
          "./result/real_data_analysis/LDL_BC.aligned",
          "./result/real_data_analysis/HDL_BC.aligned",
          "./result/real_data_analysis/BMI_BC.aligned",
          "./result/real_data_analysis/LDL_BCsub.aligned",
          "./result/real_data_analysis/HDL_BCsub.aligned",
          "./result/real_data_analysis/BMI_BCsub.aligned")
temp = 1
for(k in 1:6){
  filename = files[k]
  out_clump_SNP = read.table(filename,header=T)
  Gamma = out_clump_SNP$Gamma
  var_Gamma = out_clump_SNP$var_Gamma
  gamma = out_clump_SNP$gamma
  var_gamma = out_clump_SNP$var_gamma
  
  str_temp = strsplit(trait_vec[k],"-")
  Exposure[k] = str_temp[[1]][1]
  Outcome[k] = str_temp[[1]][2]
  #out_clump_SNP_temp = out_clump_SNP[pdx,]
 
  MR_weight_temp = MRWeight(Gamma,var_Gamma,
                            gamma,var_gamma)
  
  MRweight_result[k] <- round(MR_weight_temp[[1]],num)
  MRweight_p_value[k] = 2*pnorm(-abs(MR_weight_temp[[1]]/MR_weight_temp[[4]]),lower.tail = T)
  
  temp = temp+1
}







#three exposure
#five subtypes
bc.subtypes = c("Luminal_A","Luminal_B",
                "Luminal_B_HER2Neg",
                "HER2_Enriched",
                "Triple_Neg")
temp = 7
  for(i in 1:3){
    for(j in 1:5){
      filename = files[6+i]
      out_clump_SNP = read.table(filename,header=T)
      colnames(out_clump_SNP)[12] = "Luminal_A_log_or_meta"
      Gamma = out_clump_SNP[,paste0(bc.subtypes[j],"_log_or_meta")]
      var_Gamma = out_clump_SNP[,paste0(bc.subtypes[j],"_sd_meta")]^2
      gamma = out_clump_SNP$beta_ex
      var_gamma = out_clump_SNP$se_ex^2
      out_clump_SNP = read.table(filename,header=T)
      str_temp = strsplit(trait_vec[temp],"-")
      Exposure[temp] = str_temp[[1]][1]
      Outcome[temp] = str_temp[[1]][2]
      #out_clump_SNP_temp = out_clump_SNP[pdx,]
      
      MR_weight_temp = MRWeight(Gamma,var_Gamma,
                                gamma,var_gamma)
      
      MRweight_result[temp] <- round(MR_weight_temp[[1]],num)
      MRweight_p_value[temp] = 2*pnorm(-abs(MR_weight_temp[[1]]/MR_weight_temp[[4]]),lower.tail = T)
      
      temp = temp+1
      
    }
  }
  
result.mr.weight = data.frame(Exposure = Exposure,
                              Outcome = Outcome,
                              CasualEffect = MRweight_result,
                              p_value = MRweight_p_value,
                              stringsAsFactors = F)
result.mr.weight = result.mr.weight %>% 
  mutate(Outcome = case_when(Outcome=="BC"~"Breast Cancer",
                             Outcome=="Luminal_A"~"Luminal A",
                             Outcome=="Luminal_B"~"Luminal B",
                             Outcome=="Luminal_B_HER2Neg"~"Luminal B HER2Neg",
                             Outcome=="HER2_Enriched"~"HER2 Enriched",
                             TRUE~as.character(Outcome)),
         ) %>% 
  mutate(Exposure = factor(Exposure,levels = c("LDL","HDL","BMI")),
         Outcome = factor(Outcome, levels = c("CAD","Breast Cancer",
                          "Luminal A",
                          "Luminal B",
                          "Luminal B HER2Neg",
                          "HER2 Enriched",
                          "TN"))) %>% 
  mutate(stars = case_when(p_value<=0.0001 ~"***",
                           p_value<=0.001 ~"**",
                           p_value<=0.05 ~"*",
                           TRUE~""))
library(corrplot)
analysis.heatmap <- 
  ggplot(data = result.mr.weight, 
         aes(x = Outcome,
              y = Exposure,
              fill = CasualEffect)) + geom_tile() + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=stars), color="black", size=5)+ 
  theme_Publication()+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.3, hjust=0.3))
print(analysis.heatmap)
png(paste0("./result/real_data_analysis/analysis_heatmap.png"),height = 5, width = 6, units = "in",res =300)
print(analysis.heatmap)

dev.off()

