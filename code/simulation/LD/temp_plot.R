library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
setwd("/Users/zhangh24/GoogleDrive/MR_MA/")
source("./code/simulation/LD/theme_publication.R")
#this setting N = 60000 p =500 pthres = 5E-08
load("./result/simulation/LD/WMR_result_chr22_LDgrid.rdata")

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
  
  
tau_vec = c(0,1E-05,1E-04,1E-03,1E-02,1E-01)  
  result$tau = factor(result$tau,
                    levels = paste0("(tau = ",tau_vec,")"))

v = 1
i = 1

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
          axis.ticks.x=element_blank())+
    theme(legend.position = "none")
  print(p1)
  #p2 for confidence interval
  p2 = ggplot(result.sub,aes(x = index,y =  em_se,fill=Method))+
    geom_bar(stat="identity",position=position_dodge())+
    # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
    #               width=.2,
    #               position=position_dodge(.9))+
    theme_Publication()+
    scale_fill_Publication()+
    #coord_cartesian(ylim=c(0.75,1))+
    #geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
    facet_grid(vars(pleo_effect),vars(cau_pro),scales = "free")+
    ylab("Emirical SE")+
    ggtitle(paste0("Emirical SE"))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme(legend.position = "none")
  # print(p2)
  png(filename = paste0("./result/simulation/LD/chr22_figures/simu_result_summary_beta_",i,"_ple_",v,"LDgrid.png"),width = 15, height = 8, res=300,units  = "in")
  grid.arrange(p1,p2,p3,
               layout_matrix=rbind(c(1, 2),
                                   c(3, 3)))
  dev.off()
  
  

  r2_vec = c(0.001,0.2,0.4,0.6)
  for(r_ind in 1:length(r2_vec)){
    result.sub = result %>% filter(i_vec==i&v_vec==1&
                                     r_vec==r_ind)
    #p3 for rmse
    p3 = ggplot(result.sub,aes(x = index, y = rmse, fill = tau))+
      geom_bar(stat="identity",position=position_dodge())+
      # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
      #               width=.2,
      #               position=position_dodge(.9))+
      facet_grid(vars(pleo_effect),vars(cau_pro),scales = "free")+
      theme_Publication()+
      scale_fill_Publication()+
      ylab("RMSE")+
      ggtitle(paste0("RMSE LD r2 = ",r2_vec[r_ind]))+  
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    
    print(p3)
    #p1 for bias
    p1 = ggplot(result.sub,aes(x = index,y = bias,fill=tau))+
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
            axis.ticks.x=element_blank())+
      theme(legend.position = "none")
    print(p1)
    #p2 for confidence interval
    p2 = ggplot(result.sub,aes(x = index,y =  em_se,fill=tau))+
      geom_bar(stat="identity",position=position_dodge())+
      # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
      #               width=.2,
      #               position=position_dodge(.9))+
      theme_Publication()+
      scale_fill_Publication()+
      #coord_cartesian(ylim=c(0.75,1))+
      #geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
      facet_grid(vars(pleo_effect),vars(cau_pro),scales = "free")+
      ylab("Emirical SE")+
      ggtitle(paste0("Emirical SE"))+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      theme(legend.position = "none")
    # print(p2)
    png(filename = paste0("./result/simulation/LD/chr22_figures/simu_result_summary_beta_",i,"_ple_",v,"r_ind",r_ind,"LDgrid.png"),width = 15, height = 8, res=300,units  = "in")
    grid.arrange(p1,p2,p3,
                 layout_matrix=rbind(c(1, 2),
                                     c(3, 3)))
    dev.off()
    
  }
  
 
  
  
  
  
  v = 1
  i = 1
  
  result.sub = result %>% filter(i_vec==i&v_vec==1&
                                   tau_vec==1e-5)
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
          axis.ticks.x=element_blank())+
    theme(legend.position = "none")
  print(p1)
  #p2 for confidence interval
  p2 = ggplot(result.sub,aes(x = index,y =  em_se,fill=Method))+
    geom_bar(stat="identity",position=position_dodge())+
    # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
    #               width=.2,
    #               position=position_dodge(.9))+
    theme_Publication()+
    scale_fill_Publication()+
    #coord_cartesian(ylim=c(0.75,1))+
    #geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
    facet_grid(vars(pleo_effect),vars(cau_pro),scales = "free")+
    ylab("Emirical SE")+
    ggtitle(paste0("Emirical SE"))+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    theme(legend.position = "none")
  # print(p2)
  png(filename = paste0("./result/simulation/LD/chr22_figures/simu_result_summary_beta_",i,"_ple_",v,"tau_1e-05_LDgrid.png"),width = 15, height = 8, res=300,units  = "in")
  grid.arrange(p1,p2,p3,
               layout_matrix=rbind(c(1, 2),
                                   c(3, 3)))
  dev.off()
  
  
  