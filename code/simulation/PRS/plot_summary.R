setwd("/Users/zhangh24/GoogleDrive/MR_MA/result/simulation/PRS")
load("summary_gwas_06.rdata")
alpha_est <- result_new[[1]]  
gamma_est <- result_new[[2]]
l <- 1
plot(alpha_est[,l],gamma_est[,l])
abline(a=0,b=0.15)
points(0,0,col="red")
model <- lm(gamma_est[,l]~alpha_est[,l])
summary(model)

summary(model)



load("data_for_plot.rdata")
library(ggplot2)
data$method <- as.factor(data$method)
ggplot(data) +
  theme_Publication()+
  geom_bar(aes(x=method, y=est), stat="identity", fill="royalblue", alpha=0.7) +
  geom_errorbar(aes(x=method, ymin=est_low, ymax=est_high, width=0.2), 
                colour="firebrick2", alpha=0.9, size=0.8) +
  coord_flip() +
  # ylim(-5, 15) +
  labs(x = 'Variables', y = expression("95% CI"), 
       title = 'MR method comparasion') +
  geom_hline(yintercept = 0.15,color="red")+
  theme(plot.title=element_text(hjust=0.5, size=30),
        plot.subtitle=element_text(size=10),
        axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0), hjust=0.5, size=15),
        axis.title.y=element_text(margin=margin(t=0, r=10, b=0, l=0), vjust=0.5, size=15),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15))

