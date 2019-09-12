data = read.csv("/Users/haoyuzhang/GoogleDrive/K99_application/nih_institute_funding.csv",header=T)
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

