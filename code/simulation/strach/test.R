data = read.csv("/Users/zhangh24/GoogleDrive/K99_application/nih_institute_funding.csv",header=T)
library(dplyr)
library(ggplot2)
data_plot = data %>% select(Fiscal.Year,
                            NIH.Institute_Center,
                            Success.Rate1)
ggplot(data_plot)+
  geom_line(aes(x=Fiscal.Year,
                y=Success.Rate1,
                group=NIH.Institute_Center,
                color=NIH.Institute_Center))+
  ylab("Susscessful Rate")+
  xlab("Fiscal Year")

data_plot = data %>% select(Fiscal.Year,
                            NIH.Institute_Center,
                            Number.of.Applications.Awarded,
                            Number.of.Applications.Reviewed,
                            Success.Rate1) %>% 
  mutate(unaward = Number.of.Applications.Reviewed-Number.of.Applications.Awarded)
ggplot(data_plot)+
  geom_line(aes(x=Fiscal.Year,
                y=Number.of.Applications.Awarded,
                group=NIH.Institute_Center,
                color=NIH.Institute_Center))+
  ylab("Number of Award")+
  xlab("Fiscal Year")


data_plot = data %>% select(Fiscal.Year,
                            NIH.Institute_Center,
                            Number.of.Applications.Awarded,
                            Number.of.Applications.Reviewed,
                            Success.Rate1) %>% 
  mutate(unaward = Number.of.Applications.Reviewed-Number.of.Applications.Awarded)
ggplot(data_plot)+
  geom_line(aes(x=Fiscal.Year,
                y=Number.of.Applications.Reviewed,
                group=NIH.Institute_Center,
                color=NIH.Institute_Center))+
  ylab("Number of Application")+
  xlab("Fiscal Year")

