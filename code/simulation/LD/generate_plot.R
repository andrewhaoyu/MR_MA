library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("/Users/zhangh24/GoogleDrive/MR_MA/")
source("./code/simulation/LD/theme_publication.R")
#this setting N = 60000 p =500 pthres = 5E-08
load("./result/simulation/LD/MR_result_chr22.rdata")
result %>% filter(v_vec==1&
                    i_vec==2&method=="MRRAPs")
result %>% filter(v_vec==1&
                    i_vec==2&method=="IVW")

# 
# result.temp = result
# load("./result/simulation/LD/WMR_result_no_LD_chr22_pthres.rdata")
# result = result %>% 
#   filter(method %in%c("WMR (3e-06)","WMR (5e-04)"))
#result = rbind(result.temp,result)
beta_vec = c(0,0.2)
pleo_vec  = c(0,0.3,0.8)
cau_vec = c(0.05,0.01,0.001)
# for(i1 in 1:3){
#   for(i2 in 1:3){
result = result %>% 
  mutate(pleo_effect = case_when(v_vec ==1 ~ paste0("No Pleiotropic Effect"),
                                 v_vec==2 ~paste0(pleo_vec[2]*100,"% SNPs have Pleiotropic effect"),
                                 v_vec==3 ~paste0(pleo_vec[3]*100,"% SNPs have Pleiotropic effect")),
         cau_pro = case_when(l_vec ==1 ~ paste0(cau_vec[1]*100,"% causal SNPs"),
                             l_vec ==2 ~ paste0(cau_vec[2]*100,"% causal SNPs"),
                             l_vec ==3 ~ paste0(cau_vec[3]*100,"% causal SNPs")),
         index = 1) %>% 
  #rename(Method = method) %>% 
  mutate(pleo_effect = factor(pleo_effect,
                                  levels = c("No Pleiotropic Effect",paste0(
                                                  pleo_vec[2:3]*100,"% SNPs have Pleiotropic effect")))) %>% 
  filter(method!="MR-Egger")


p3 = ggplot(result,aes(x = index, y = rmse, fill = method))+
  geom_bar(stat="identity",position=position_dodge())+
  # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
  #               width=.2,
  #               position=position_dodge(.9))+
  facet_grid(vars(cau_pro),vars(pleo_effect))+
  theme_Publication()+
  scale_fill_Publication()+
  ylab("RMSE")+
  ggtitle(paste0("RMSE"))+  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p1 = ggplot(result,aes(x = index,y = bias,fill=Method))+
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
p2 = ggplot(result,aes(x = index,y =  cover,fill=Method))+
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
load("./result/simulation/LD/WMR_result_no_LD_chr22_pthres.rdata")
pthres = c(3E-06,1E-05,5E-05,1E-04,5E-04,1E-03)       
result = result %>% 
  mutate(beta_effect = case_when(i_vec ==1 ~ paste0(beta_vec[1]*100,"% h2 of Y mediated by M"),
                                 i_vec==2 ~paste0(beta_vec[2]^2*100,"% h2 of Y mediated by M"),
                                 i_vec==3 ~paste0(beta_vec[3]^2*100,"% h2 of Y mediated by M"))) %>%
  mutate(index = 1,
         beta_effect = factor(beta_effect,
                              levels = paste0(beta_vec^2*100,"% h2 of Y mediated by M")),
         method = factor(method, levels = paste0("WMR (",pthres[length(pthres):1], ")"))
  ) %>% 
  rename(Method=method)

p3 = ggplot(result,aes(x = index, y = rmse, fill = Method))+
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
p1 = ggplot(result,aes(x = index,y = bias,fill=Method))+
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


png(filename = paste0("./result/simulation/LD/chr22_wmrsimu_noLD_result_rmse.png"),width = 12, height = 8, res=300,units  = "in")
print(p3)
dev.off()

png(filename = paste0("./result/simulation/LD/chr22_wmrsimu_noLD_result_bias.png"),width = 12, height = 8, res=300,units  = "in")
print(p1)
dev.off()




#wmr result under different p_thres
load("./result/simulation/LD/WMR_result_chr22_pthres.rdata")
pthres = c(3E-06,1E-05,5E-05,1E-04,5E-04,1E-03)       
result = result %>% 
  mutate(beta_effect = case_when(i_vec ==1 ~ paste0(beta_vec[1]*100,"% h2 of Y mediated by M"),
                                 i_vec==2 ~paste0(beta_vec[2]^2*100,"% h2 of Y mediated by M"),
                                 i_vec==3 ~paste0(beta_vec[3]^2*100,"% h2 of Y mediated by M"))) %>%
  mutate(index = 1,
         beta_effect = factor(beta_effect,
                              levels = paste0(beta_vec^2*100,"% h2 of Y mediated by M")),
         method = factor(method, levels = paste0("WMR (",pthres[length(pthres):1], ")"))
  ) %>% 
  rename(Method=method)

p3 = ggplot(result,aes(x = index, y = rmse, fill = Method))+
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
p1 = ggplot(result,aes(x = index,y = bias,fill=Method))+
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


png(filename = paste0("./result/simulation/LD/chr22_wmrsimu_strongLD_result_rmse.png"),width = 12, height = 8, res=300,units  = "in")
print(p3)
dev.off()

png(filename = paste0("./result/simulation/LD/chr22_wmrsimu_strongLD_result_bias.png"),width = 12, height = 8, res=300,units  = "in")
print(p1)
dev.off()








load("./result/simulation/LD/WMR_result_chr22_pthres.rdata")
result = result %>% filter(method =="WMR (1e-06)") %>% 
  mutate(method = paste0("Strong LD ",method))
result.temp = result
load("./result/simulation/LD/WMR_result_no_LD_chr22_pthres.rdata")
pthres = c(3E-06,1E-05,5E-05,1E-04,5E-04,1E-03)       
result = result %>% 
  mutate(beta_effect = case_when(i_vec ==1 ~ paste0(beta_vec[1]*100,"% h2 of Y mediated by M"),
                                 i_vec==2 ~paste0(beta_vec[2]^2*100,"% h2 of Y mediated by M"),
                                 i_vec==3 ~paste0(beta_vec[3]^2*100,"% h2 of Y mediated by M"))) %>%
  mutate(index = 1,
         beta_effect = factor(beta_effect,
                              levels = paste0(beta_vec^2*100,"% h2 of Y mediated by M")),
         method = factor(method, levels = paste0("WMR (",pthres[length(pthres):1], ")"))
  ) %>% 
  rename(Method=method)

p3 = ggplot(result,aes(x = index, y = rmse, fill = Method))+
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
p1 = ggplot(result,aes(x = index,y = bias,fill=Method))+
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


png(filename = paste0("./result/simulation/LD/chr22_wmrsimu_strongLD_result_rmse.png"),width = 12, height = 8, res=300,units  = "in")
print(p3)
dev.off()

png(filename = paste0("./result/simulation/LD/chr22_wmrsimu_strongLD_result_bias.png"),width = 12, height = 8, res=300,units  = "in")
print(p1)
dev.off()
