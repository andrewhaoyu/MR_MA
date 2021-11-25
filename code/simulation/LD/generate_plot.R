library(dplyr)
library(ggplot2)
library(ggpubr)

setwd("/Users/zhangh24/GoogleDrive/MR_MA/")
source("./code/simulation/LD/theme_publication.R")
#this setting N = 60000 p =500 pthres = 5E-08
load("./result/simulation/LD/MR_result_chr22.rdata")
result.stan = result
load("./result/simulation/LD/wmr_simu_result_noLD_pthres.rdata")
result.wmr.nold = result
result.wmr.nold.sub = result.wmr.nold %>% 
  filter(method%in%paste0("WMR (p<",c(5e-08,1e-04),")"))
result = rbind(result.stan,result.wmr.nold.sub)

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
  mutate(pleo_effect = case_when(v_vec ==1 ~ paste0("No pleiotropic effect"),
                                 v_vec==2 ~paste0("Mild pleiotropic effect"),
                                 v_vec==3 ~paste0("High pleiotropic effect")),
         cau_pro = case_when(l_vec ==1 ~ paste0(cau_vec[1]*100,"% causal SNPs"),
                             l_vec ==2 ~ paste0(cau_vec[2]*100,"% causal SNPs"),
                             l_vec ==3 ~ paste0(cau_vec[3]*100,"% causal SNPs")),
         index = 1) %>% 
  #rename(Method = method) %>% 
  mutate(pleo_effect = factor(pleo_effect,
                                  levels = c("No pleiotropic effect","Mild pleiotropic effect","High pleiotropic effect"))) %>% 
  filter(method!="MR-Egger") %>% 
  rename(Method = method)
for(i in 1:2){
  result.sub = result %>% filter(i_vec==i&v_vec==1)
  #p3 for rmse
  p3 = ggplot(result.sub,aes(x = index, y = rmse, fill = Method))+
    geom_bar(stat="identity",position=position_dodge())+
    # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
    #               width=.2,
    #               position=position_dodge(.9))+
    facet_grid(vars(pleo_effect),vars(cau_pro),scales = "free")+
    theme_Publication()+
    scale_fill_Publication()+
    ylab("RMSE")+
    ggtitle(paste0("RMSE"))+  
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  print(p3)
  #p1 for bias
  p1 = ggplot(result.sub,aes(x = index,y = bias,fill=Method))+
    geom_bar(stat="identity",position=position_dodge())+
    theme_Publication()+
    scale_fill_Publication()+
    #ylim(c(-0.1,0.1))+
    ylab("Bias estimate")+
    xlab(NULL)+
    facet_grid(vars(pleo_effect),vars(cau_pro),scales = "free")+
    ggtitle(paste0("Bias"))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  print(p1)
  #p2 for confidence interval
  p2 = ggplot(result.sub,aes(x = index,y =  cover,fill=Method))+
    geom_bar(stat="identity",position=position_dodge())+
    # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
    #               width=.2,
    #               position=position_dodge(.9))+
    theme_Publication()+
    scale_fill_Publication()+
    coord_cartesian(ylim=c(0.75,1))+
    geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
    facet_grid(vars(pleo_effect),vars(cau_pro),scales = "free")+
    ylab("Coverage")+
    ggtitle(paste0("Coverage"))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  print(p2)
  pthres = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02)
  result.wmr.nold.plot= result.wmr.nold %>% 
    mutate(pleo_effect = case_when(v_vec ==1 ~ paste0("No pleiotropic effect"),
                                   v_vec==2 ~paste0("Mild pleiotropic effect"),
                                   v_vec==3 ~paste0("High pleiotropic effect")),
           cau_pro = case_when(l_vec ==1 ~ paste0(cau_vec[1]*100,"% causal SNPs"),
                               l_vec ==2 ~ paste0(cau_vec[2]*100,"% causal SNPs"),
                               l_vec ==3 ~ paste0(cau_vec[3]*100,"% causal SNPs")),
           index = 1) %>% 
    #rename(Method = method) %>% 
    mutate(pleo_effect = factor(pleo_effect,
                                levels = c("No pleiotropic effect","Mild pleiotropic effect","High pleiotropic effect")),
           method = factor(method,levels = paste0("WMR (p<",pthres[c(length(pthres):1)],")"))) %>% 
    rename(Method = method)
  result.sub = result.wmr.nold.plot %>% 
    filter(i_vec==i&v_vec==1)
  
  p4 = ggplot(result.sub,aes(x = index,y = rmse,fill=Method))+
    geom_bar(stat="identity",position=position_dodge())+
    theme_Publication()+
    scale_fill_Publication()+
    #ylim(c(-0.1,0.1))+
    ylab("RMSE estimate")+
    xlab(NULL)+
    facet_grid(vars(pleo_effect),vars(cau_pro),scales = "free")+
    ggtitle(paste0("RMSE"))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  print(p4)
  png(filename = paste0("./result/simulation/LD/chr22_figures/simu_result_rmse_beta_",i,"_ple_",v,".png"),width = 12, height = 8, res=300,units  = "in")
  print(p3)
  dev.off()
  
  png(filename = paste0("./result/simulation/LD/chr22_figures/chr22_simu_result_bias_beta_",i,"_ple_",v,".png"),width = 12, height = 8, res=300,units  = "in")
  print(p1)
  dev.off()
  png(filename = paste0("./result/simulation/LD/chr22_figures/chr22_simu_result_coverage_beta_",i,"_ple_",v,".png"),width = 12, height = 8, res=300,units  = "in")
  print(p2)
  dev.off()
  png(filename = paste0("./result/simulation/LD/chr22_figures/chr22_simu_resul_wmr_beta_",i,"_ple_",v,".png"),width = 12, height = 8, res=300,units  = "in")
  print(p4)
  dev.off()
  
}

#compare low LD and high LD settings
load("./result/simulation/LD/wmr_simu_result_strongLD_pthres.rdata")

result.wmr.sld = result
result.wmr.sld.sub = result.wmr.sld %>% 
  filter(method%in%paste0("WMR (p<",c(5e-08,1e-04),")")) %>% 
  mutate(method=paste0("High LD ",method))
load("./result/simulation/LD/wmr_simu_result_noLD_pthres.rdata")
result.wmr.nold = result
result.wmr.nold.sub = result.wmr.nold %>% 
  filter(method%in%paste0("WMR (p<",c(5e-08,1e-04),")")) %>% 
  mutate(method=paste0("Low LD ",method))
result = rbind(result.wmr.nold.sub,result.wmr.sld.sub)
result = result %>% 
  mutate(pleo_effect = case_when(v_vec ==1 ~ paste0("No pleiotropic effect"),
                                 v_vec==2 ~paste0("Mild pleiotropic effect"),
                                 v_vec==3 ~paste0("High pleiotropic effect")),
         cau_pro = case_when(l_vec ==1 ~ paste0(cau_vec[1]*100,"% causal SNPs"),
                             l_vec ==2 ~ paste0(cau_vec[2]*100,"% causal SNPs"),
                             l_vec ==3 ~ paste0(cau_vec[3]*100,"% causal SNPs")),
         index = 1) %>% 
  #rename(Method = method) %>% 
  mutate(pleo_effect = factor(pleo_effect,
                              levels = c("No pleiotropic effect","Mild pleiotropic effect","High pleiotropic effect"))) %>% 
  filter(method!="MR-Egger") %>% 
  rename(Method = method)
result.sub = result %>% filter(i_vec==1&v_vec==1)
p3 = ggplot(result.sub,aes(x = index, y = rmse, fill = Method))+
  geom_bar(stat="identity",position=position_dodge())+
  # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
  #               width=.2,
  #               position=position_dodge(.9))+
  facet_grid(vars(pleo_effect),vars(cau_pro),scales = "free")+
  theme_Publication()+
  scale_fill_Publication()+
  ylab("RMSE")+
  ggtitle(paste0("RMSE"))+  
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

print(p3)
png(filename = paste0("./result/simulation/LD/chr22_figures/LDcomparasion_rmse_beta_",i,"_ple_",v,".png"),width = 12, height = 8, res=300,units  = "in")
print(p3)
dev.off()
