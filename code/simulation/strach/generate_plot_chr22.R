library(tidyverse)
setwd("/Users/zhangh24/GoogleDrive/MR_MA/result/simulation")
data = read.csv("simulatation_22_temp.csv")
data_plot = data %>% 
  filter(Pleiotropic==1) %>% 
  mutate(Index = 1) %>% 
  mutate(Beta = as.numeric(Beta))

source("/Users/zhangh24/GoogleDrive/MR_MA/code/simulation/LD/theme_publication.R")

p1 = ggplot(data_plot,aes(Index,y = Beta,fill=Method)) + 
  geom_bar(stat="identity",
           position = position_dodge()) + 
  facet_grid(cols=vars(Cau))+
  theme_Publication()+
  geom_hline(yintercept = 0.15, linetype = "dashed")+
  scale_fill_Publication()+
  xlab(NULL)+
 # scale_y_continuous(limits = c(0.10, 0.20))+
  theme(axis.text.x=element_blank())
  
p2 = ggplot(data_plot,aes(Index,y = SE,fill=Method)) + 
  geom_bar(stat="identity",
           position = position_dodge()) + 
  facet_grid(cols=vars(Cau))+
  theme_Publication()+
  scale_fill_Publication()+
  xlab(NULL)+
  theme(axis.text.x=element_blank())


data_plot = data %>% 
  filter(Pleiotropic==2) %>% 
  mutate(Index = 1) %>% 
  mutate(Beta = as.numeric(Beta))

p3 = ggplot(data_plot,aes(Index,y = Beta,fill=Method)) + 
  geom_bar(stat="identity",
           position = position_dodge()) + 
  facet_grid(cols=vars(Cau))+
  theme_Publication()+
  geom_hline(yintercept = 0.15, linetype = "dashed")+
  scale_fill_Publication()+
  xlab(NULL)+
  #scale_y_continuous(limits = c(0.10, 0.20))+
  theme(axis.text.x=element_blank())

p4 = ggplot(data_plot,aes(Index,y = SE,fill=Method)) + 
  geom_bar(stat="identity",
           position = position_dodge()) + 
  facet_grid(cols=vars(Cau))+
  theme_Publication()+
  scale_fill_Publication()+
  xlab(NULL)+
  theme(axis.text.x=element_blank())

png(filename = "beta_no_pleo",width=10,height = 8,units = "in",res=300)
print(p1)
dev.off()

png(filename = "se_no_pleo",width=10,height = 8,units = "in",res=300)
print(p2)
dev.off()

png(filename = "beta_with_pleo",width=10,height = 8,units = "in",res=300)
print(p3)
dev.off()

png(filename = "se_with_pleo",width=10,height = 8,units = "in",res=300)
print(p4)
dev.off()
