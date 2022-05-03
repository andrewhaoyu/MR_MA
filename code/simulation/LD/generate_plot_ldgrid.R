library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
setwd("/Users/zhangh24/GoogleDrive/MR_MA/")
source("./code/simulation/LD/theme_publication.R")
#this setting N = 60000 p =500 pthres = 5E-08
load("./result/simulation/LD/WMR_result_chr22_grid.rdata")
head(result)
#subset result to WMR
result = result %>% filter(method=="WMR")
#v_vec = 1 as no pleiotropic effect
result$v_vec = 1
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

r2_vec = c(0.001,0.2,0.4,0.6,0.8,1)
pthres = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02,1E-01,1)
#fixed r_vec as 0.4 (r_vec=3) and see the effect of p-value
for(r in c(0.001,0.4)){
  for(i in 1:2){
    for(l in 1:3){
      result.sub = result %>% filter(i_vec==i&r_vec==r&
                                       l_vec ==l) %>% 
        mutate(p_vec = factor(p_vec,
                              levels = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02,1E-01,1)))
      max.y = max(c(result.sub$rmse,result.sub$em_se))
      p3 = ggplot(result.sub,aes(x = index, y = rmse, fill = p_vec))+
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
              axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="P-value"))
      
      #print(p3)
      #p1 for bias
      p1 = ggplot(result.sub,aes(x = index,y = bias,fill=p_vec))+
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
        theme(legend.position = "none")+
        scale_y_continuous(limits = c(NA,max.y))+
        guides(fill=guide_legend(title="P-value"))
      #print(p1)
      #p2 for confidence interval
      p2 = ggplot(result.sub,aes(x = index,y =  em_se,fill=p_vec))+
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
        theme(legend.position = "none")+
        scale_y_continuous(limits = c(NA,max.y))+
        guides(fill=guide_legend(title="P-value"))
      # print(p2)
      png(filename = paste0("./result/simulation/LD/chr22_figures/fidxed_r2_different_p_value_r_",r,"_beta_",i,"_cau_",l,"_LDgrid.png"),width = 15, height = 8, res=300,units  = "in")
      grid.arrange(p1,p2,p3,
                   layout_matrix=rbind(c(1, 2),
                                       c(3, 3)),
                   top=textGrob(paste0("Beta = ",beta_vec[i],", clumping R2 = ",r),gp=gpar(fontsize=24,font = 2)))
      dev.off()
      
  }
      
    }
  
}

r2_vec = c(0.001,0.2,0.4,0.6,0.8,1)
pthres = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02,1E-01,1)
p = 5
#fixed p-value as 1E-04 and see the effect of R2
for(i in 1:2){
    for(l in 1:3){
      result.sub = result %>% filter(i_vec==i&p_vec==pthres[p]&
                                       l_vec ==l) %>% 
        mutate(r_vec = factor(r_vec,
                              levels = c(0.001,0.2,0.4,0.6,0.8,1)))
      max.y = max(c(result.sub$rmse,result.sub$em_se))
      p3 = ggplot(result.sub,aes(x = index, y = rmse, fill = r_vec))+
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
              axis.ticks.x=element_blank())+
        guides(fill=guide_legend(title="Clumping-R2"))
      
      #print(p3)
      #p1 for bias
      p1 = ggplot(result.sub,aes(x = index,y = bias,fill=r_vec))+
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
        theme(legend.position = "none")+
        scale_y_continuous(limits = c(NA,max.y))+
        guides(fill=guide_legend(title="Clumping-R2"))
      #print(p1)
      #p2 for confidence interval
      p2 = ggplot(result.sub,aes(x = index,y =  em_se,fill=r_vec))+
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
        theme(legend.position = "none")+
        scale_y_continuous(limits = c(NA,max.y))+
        guides(fill=guide_legend(title="Clumping-R2"))
      # print(p2)
      png(filename = paste0("./result/simulation/LD/chr22_figures/fidxed_p_value_different_clumpingr2_beta_",i,"_cau_",l,"_LDgrid.png"),width = 15, height = 8, res=300,units  = "in")
      grid.arrange(p1,p2,p3,
                   layout_matrix=rbind(c(1, 2),
                                       c(3, 3)),
                   top=textGrob(paste0("Beta = ",beta_vec[i]," p-value = ",pthres[p]),gp=gpar(fontsize=24,font = 2)))
      dev.off()
      
    }
    
  }
  



#check the performance of other approaches
load("./result/simulation/LD/WMR_result_chr22_grid.rdata")
head(result)
#subset result to WMR
#v_vec = 1 as no pleiotropic effect
result$v_vec = 1
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

r2_vec = c(0.001)
pthres = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02,1E-01,1)
method_vec = c("Raps","divw")
for(m in 1:length(method_vec)){
  for(r in c(0.001)){
    for(i in 1:2){
      for(l in 1:3){
        method = method_vec[m]
        result.sub = result %>% filter(i_vec==i&r_vec==r&
                                         l_vec ==l&
                                         Method==method) %>% 
          mutate(p_vec = factor(p_vec,
                                levels = c(5E-08,1E-07,1E-06,1E-05,1E-04,1E-03,1E-02,1E-01,1)))
        max.y = max(c(result.sub$rmse,result.sub$em_se))
        p3 = ggplot(result.sub,aes(x = index, y = rmse, fill = p_vec))+
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
                axis.ticks.x=element_blank())+
          guides(fill=guide_legend(title="P-value"))
        
        #print(p3)
        #p1 for bias
        p1 = ggplot(result.sub,aes(x = index,y = bias,fill=p_vec))+
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
          theme(legend.position = "none")+
          scale_y_continuous(limits = c(NA,max.y))+
          guides(fill=guide_legend(title="P-value"))
        #print(p1)
        #p2 for confidence interval
        p2 = ggplot(result.sub,aes(x = index,y =  em_se,fill=p_vec))+
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
          theme(legend.position = "none")+
          scale_y_continuous(limits = c(NA,max.y))+
          guides(fill=guide_legend(title="P-value"))
        # print(p2)
        png(filename = paste0("./result/simulation/LD/chr22_figures/",method,",_r2_different_p_value_r_",r,"_beta_",i,"_cau_",l,"_LDgrid.png"),width = 15, height = 8, res=300,units  = "in")
        grid.arrange(p1,p2,p3,
                     layout_matrix=rbind(c(1, 2),
                                         c(3, 3)),
                     top=textGrob(paste0(method," Beta = ",beta_vec[i],", clumping R2 = ",r),gp=gpar(fontsize=24,font = 2)))
        dev.off()
        
      }
      
    }
    
  }
}
