setwd("/Users/zhangh24/GoogleDrive/MR_MA/")
source("./code/simulation/LD/theme_publication.R")
#this setting N = 60000 p =500 pthres = 5E-08
load("./result/simulation/LD/MR_result_chr22.rdata")
result.temp = result
load("./result/simulation/LD/WMR_result_chr22_i14.rdata")
result = rbind(result.temp,result)
beta_vec = round(c(1,0.5,0),2)
pleo_vec  = c(1,0.5,0.25)
library(dplyr)
library(ggplot2)
library(ggpubr)
# for(i1 in 1:3){
#   for(i2 in 1:3){
result = result %>% 
  mutate(beta_effect = case_when(i_vec ==1 ~ paste0(beta_vec[1]*100,"% h2 of Y mediated by M"),
                                 i_vec==2 ~paste0(beta_vec[2]^2*100,"% h2 of Y mediated by M"),
                                 i_vec==3 ~paste0(beta_vec[3]^2*100,"% h2 of Y mediated by M"))) %>%
  mutate(index = 1,
  beta_effect = factor(beta_effect,
                              levels = paste0(beta_vec^2*100,"% h2 of Y mediated by M"),
                              ),
  )
  

p3 = ggplot(result,aes(x = index, y = rmse, fill = method))+
  geom_bar(stat="identity",position=position_dodge())+
  # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
  #               width=.2,
  #               position=position_dodge(.9))+
  facet_wrap(vars(beta_effect))+
               theme_Publication()+
               scale_fill_Publication()+
               ylab("RMSE")+
               ggtitle(paste0("RMSE"))+  
  theme(axis.title.x=element_blank(),
                                               axis.text.x=element_blank(),
                                               axis.ticks.x=element_blank())
             p1 = ggplot(result,aes(x = index,y = bias,fill=method))+
               geom_bar(stat="identity",position=position_dodge())+
               theme_Publication()+
               scale_fill_Publication()+
               #ylim(c(-0.1,0.1))+
               ylab("Bias estimate")+
               xlab(NULL)+
               facet_wrap(vars(beta_effect))+
               ggtitle(paste0("Bias"))+
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
             p2 = ggplot(result,aes(x = index,y =  cover,fill=method))+
               geom_bar(stat="identity",position=position_dodge())+
               # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
               #               width=.2,
               #               position=position_dodge(.9))+
               theme_Publication()+
               scale_fill_Publication()+
               coord_cartesian(ylim=c(0.75,1))+
               geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
               facet_wrap(vars(beta_effect))+
               ylab("Coverage")+
               ggtitle(paste0("Coverage"))+
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
             png(filename = paste0("./result/simulation/LD/chr22_simu_result_rmse.png"),width = 12, height = 8, res=300,units  = "in")
             print(p3)
             dev.off()
             
             png(filename = paste0("./result/simulation/LD/chr22_simu_result_bias.png"),width = 12, height = 8, res=300,units  = "in")
             print(p1)
             dev.off()
             
             # png(filename = paste0("./result/simulation/LD/test_simu_result_coverage.png"),width = 12, height = 8, res=300,units  = "in")
             # print(p2)
             # dev.off()

             
#wmr result under different p_thres
load("./result/simulation/LD/wmr_simu_result_com_np_pthres.rdata")
             
result = result %>% 
  mutate(beta_effect = case_when(i1_vec ==1 ~ paste0(beta_vec[1]*100,"% h2 of Y mediated by M"),
                                 i1_vec==2 ~paste0(beta_vec[2]^2*100,"% h2 of Y mediated by M"),
                                 i1_vec==3 ~paste0(beta_vec[3]^2*100,"% h2 of Y mediated by M"))) %>% 
  mutate(cau_pro = case_when(i2_vec ==1 ~ paste0(pleo_vec[1]*100,"% causal SNPs shared"),
                             i2_vec==2 ~paste0(pleo_vec[2]*100,"% causal SNPs shared"),
                             i2_vec==3 ~paste0(pleo_vec[3]*100,"% causal SNPs shared"))) %>% 
  mutate(index = 1,
          beta_effect = factor(beta_effect,
          levels = paste0(beta_vec^2*100,"% h2 of Y mediated by M"),
          ),
          cau_pro = factor(cau_pro,
          levels = paste0(pleo_vec*100,"% causal SNPs shared")
          ))
result.sub= result %>% 
  filter(i1_vec==1&i2_vec==1)
pthres = c(1,0.1,0.01,1E-03,1E-4,1E-05)
p1 = ggplot(result,aes(x = -log10(method),y = bias))+
  geom_bar(stat="identity")+
  theme_Publication()+
  #scale_fill_Publication()+
  #ylim(c(-0.1,0.1))+
  ylab("Bias estimate")+
  xlab("-log10(Pvalue)")+
  facet_grid(vars(cau_pro),vars(beta_effect))+
  ggtitle(paste0("Bias"))+
  scale_x_continuous(breaks = -log10(pthres))
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())


p3 = ggplot(result,aes(x = -log10(method),y = rmse))+
  geom_bar(stat="identity",position=position_dodge())+
  theme_Publication()+
  #scale_fill_Publication()+
  #ylim(c(-0.1,0.1))+
  ylab("RMSE estimate")+
  xlab("-log10(Pvalue)")+
  facet_grid(vars(cau_pro),vars(beta_effect))+
  ggtitle(paste0("RMSE"))+
  scale_x_continuous(breaks = -log10(pthres))

png(filename = paste0("./result/simulation/LD/test_simu_result_np_pthres_bias.png"),width = 12, height = 8, res=300,units  = "in")
print(p1)
dev.off()
png(filename = paste0("./result/simulation/LD/test_simu_result_np_pthres_rmse.png"),width = 12, height = 8, res=300,units  = "in")
print(p3)
dev.off()
  result.sub = result %>% filter(method ==1E-05)  %>% 
    mutate(method = "WMR (1E-05)")
#this setting N = 2000 p =1000 pthres for IVW = 1E-05
load("./result/simulation/LD/wmr_simu_result_com_np.rdata")
beta_vec = round(c(1,0.5,0),2)
pleo_vec  = c(1,0.5,0.25)
library(dplyr)
library(ggplot2)
library(ggpubr)
             # for(i1 in 1:3){
             #   for(i2 in 1:3){
result = result %>% 
  mutate(beta_effect = case_when(i1_vec ==1 ~ paste0(beta_vec[1]*100,"% h2 of Y mediated by M"),
                                 i1_vec==2 ~paste0(beta_vec[2]^2*100,"% h2 of Y mediated by M"),
                                 i1_vec==3 ~paste0(beta_vec[3]^2*100,"% h2 of Y mediated by M"))) %>% 
  mutate(cau_pro = case_when(i2_vec ==1 ~ paste0(pleo_vec[1]*100,"% causal SNPs shared"),
                             i2_vec==2 ~paste0(pleo_vec[2]*100,"% causal SNPs shared"),
                             i2_vec==3 ~paste0(pleo_vec[3]*100,"% causal SNPs shared"))) %>% 
  mutate(index = 1,
         beta_effect = factor(beta_effect,
                              levels = paste0(beta_vec^2*100,"% h2 of Y mediated by M"),
         ),
         cau_pro = factor(cau_pro,
                          levels = paste0(pleo_vec*100,"% causal SNPs shared")
         ))
result = rbind(result,result.sub)

p3 = ggplot(result,aes(x = index, y = rmse, fill = method))+
               geom_bar(stat="identity",position=position_dodge())+
               # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
               #               width=.2,
               #               position=position_dodge(.9))+
               facet_grid(vars(cau_pro),vars(beta_effect))+
               theme_Publication()+
               scale_fill_Publication()+
               ylab("RMSE")+
               ggtitle(paste0("RMSE"))+  
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
             p1 = ggplot(result,aes(x = index,y = bias,fill=method))+
               geom_bar(stat="identity",position=position_dodge())+
               theme_Publication()+
               scale_fill_Publication()+
               #ylim(c(-0.1,0.1))+
               ylab("Bias estimate")+
               xlab(NULL)+
               facet_grid(vars(cau_pro),vars(beta_effect))+
               ggtitle(paste0("Bias"))+
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
             p2 = ggplot(result,aes(x = index,y =  cover,fill=method))+
               geom_bar(stat="identity",position=position_dodge())+
               # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
               #               width=.2,
               #               position=position_dodge(.9))+
               theme_Publication()+
               scale_fill_Publication()+
               coord_cartesian(ylim=c(0.75,1))+
               geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
               facet_grid(vars(cau_pro),vars(beta_effect))+
               ylab("Coverage")+
               ggtitle(paste0("Coverage"))+
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank())
             png(filename = paste0("./result/simulation/LD/test_simu_result_np_rmse.png"),width = 12, height = 8, res=300,units  = "in")
             print(p3)
             dev.off()
             
             png(filename = paste0("./result/simulation/LD/test_simu_result_np_bias.png"),width = 12, height = 8, res=300,units  = "in")
             print(p1)
             dev.off()
             
             png(filename = paste0("./result/simulation/LD/test_simu_result_np_coverage.png"),width = 12, height = 8, res=300,units  = "in")
             print(p2)
             dev.off()
             
             
             
             