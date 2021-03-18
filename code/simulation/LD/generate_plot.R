setwd("/Users/zhangh24/GoogleDrive/MR_MA/")
source("./code/simulation/LD/theme_publication.R")
load("./result/simulation/LD/MR_result_summary.rdata")
library(dplyr)
cau_vec = c(0.05,0.01,0.001)

for(l in 2:3){
  MR_result_temp= MR_result_com %>% 
    filter(l_vec==l) %>% 
    mutate(method = factor(method,levels = c("IVW","MR-Egger","MR-Median","MR-PRESSO","Raps","MR-Weight"))) %>% 
    group_by(method) %>% 
    summarize(mean_est = mean(est,na.rm =T),
              sd_est = sd(est,na.rm=T),
              mean_cover = mean(cover,na.rm =T),
              nonna_count = sum(!is.na(est))) %>% 
    mutate(sample_size = "60000") 
  
  library(ggplot2)
  
  p1 = ggplot(MR_result_temp,aes(x = sample_size,y = mean_est,fill=method))+
    geom_bar(stat="identity",color="black",position=position_dodge())+
    geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
                  width=.2,
                  position=position_dodge(.9))+
    theme_Publication()+
    scale_fill_Publication()+
    geom_hline(yintercept = 0.15, col="red",linetype ="dashed")+
    ylab("Beta estimate")+
    ggtitle(paste0("Causal SNPs Proportion = ",cau_vec[l]))+
    theme(legend.position = "none")
  
  p2 = ggplot(MR_result_temp,aes(x = sample_size,y =  mean_cover,fill=method))+
    geom_bar(stat="identity",color="black",position=position_dodge())+
    # geom_errorbar(aes(ymin=mean_est-sd_est,ymax=mean_est+sd_est),
    #               width=.2,
    #               position=position_dodge(.9))+
    theme_Publication()+
    scale_fill_Publication()+
    coord_cartesian(ylim=c(0.75,1))+
    geom_hline(yintercept = 0.95, col="red",linetype ="dashed")+
    ylab("Coverage")+
    ggtitle(paste0("Causal SNPs Proportion = ",cau_vec[l]))
  
  library(ggpubr)
  
  p = ggarrange(p1,p2)
  png(filename = paste0("MR_result_rho_",l,".png"),width = 12, height = 8, res=300,units  = "in")
  print(p)
  dev.off()
  
  
}
