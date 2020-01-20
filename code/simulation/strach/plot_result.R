load("./result/simulation/data_plot_fixed_g_0.01_result.rdata") 
library(moments)
data1 <- as.data.frame(data[[1]])
library(dplyr)
data1 = data1 %>% mutate(
  TSN = dnorm(TwoStage_est,
               mean=mean(TwoStage_est),
               sd=sd(TwoStage_est)),
  IVWN = dnorm(IVW_est,
               mean=mean(IVW_est),
               sd=sd(IVW_est)),
  IVWmeta_N = dnorm(IVW_meta_est,
               mean=mean(IVW_meta_est),
               sd=sd(IVW_meta_est)),
  IVWs_N = dnorm(IVWs_est,
                    mean=mean(IVWs_est),
                    sd=sd(IVWs_est))
)

skewness

kurtosis

ske_vec <- rep(0,4)
ske_vec_n <- rep(0,4)
kur_vec <- rep(0,4)
kur_vec_n <- rep(0,4)
for(k in 1:4){
  ske_vec[k] <- skewness(data1[,k])
  kur_vec[k] <- kurtosis(data1[,k])
  temp_n <- rnorm(10000,
                  mean = mean(data1[,k]),
                  sd = sd(data1[,k]))
  ske_vec_n[k] <- skewness(temp_n)
  kur_vec_n[k] <- kurtosis(temp_n)
  }
result <- cbind(ske_vec,kur_vec,
                ske_vec_n,kur_vec_n)

row.names(result) <- c("Two Stage",
                       "IVW","IVW meta",
                       "IVW summary")
temp_n <- rnorm(100000,
                0,
                10)
temp_t <- rt(1000000,
             df=5.363636)
kurtosis(temp_n)
kurtosis(temp_t)


library(ggplot2)
p1 = ggplot(data1)+
  geom_histogram(aes(x=TwoStage_est,y=..density..))+
  geom_line(aes(x=TwoStage_est,y=TSN),color="red")+
  theme_Publication()+
  xlab("Two Stage Estimate")+
  ylim(limits=c(0,1.8))+
  xlim(limits=c(-2,2))

p2 = ggplot(data1)+
  geom_histogram(aes(x=IVW_est,y=..density..))+
  geom_line(aes(x=IVW_est,y=IVWN),color="red")+
  theme_Publication()+
  xlab("IVW")+
  ylim(limits=c(0,1.8))+
  xlim(limits=c(-2,2))


p3 = ggplot(data1)+
  geom_histogram(aes(x=IVW_meta_est,y=..density..))+
  geom_line(aes(x=IVW_meta_est,y=IVWmeta_N),color="red")+
  theme_Publication()+
  xlab("IVW meta Estimate")+
  ylim(limits=c(0,1.8))+
  xlim(limits=c(-2,2))

p4 = ggplot(data1)+
  geom_histogram(aes(x=IVWs_est,y=..density..))+
  geom_line(aes(x=IVWs_est,y=IVWs_N),color="red")+
  theme_Publication()+
  xlab("IVW summary")+
  ylim(limits=c(0,1.8))+
  xlim(limits=c(-2,2))
library(gridExtra)
grid.arrange(p1,p2,p4,nrow=2)











data1 <- as.data.frame(data[[2]])
library(dplyr)
data1 = data1 %>% mutate(
  TSN = dnorm(TwoStage_est,
              mean=mean(TwoStage_est),
              sd=sd(TwoStage_est)),
  IVWN = dnorm(IVW_est,
               mean=mean(IVW_est),
               sd=sd(IVW_est)),
  IVWmeta_N = dnorm(IVW_meta_est,
                    mean=mean(IVW_meta_est),
                    sd=sd(IVW_meta_est)),
  IVWs_N = dnorm(IVWs_est,
                 mean=mean(IVWs_est),
                 sd=sd(IVWs_est))
)

p1 = ggplot(data1)+
  geom_histogram(aes(x=TwoStage_est,y=..density..))+
  geom_line(aes(x=TwoStage_est,y=TSN),color="red")+
  theme_Publication()+
  xlab("Two Stage Estimate")+
  ylim(limits=c(0,4))

p2 = ggplot(data1)+
  geom_histogram(aes(x=IVW_est,y=..density..))+
  geom_line(aes(x=IVW_est,y=IVWN),color="red")+
  theme_Publication()+
  xlab("IVW")+
  ylim(limits=c(0,4))


p3 = ggplot(data1)+
  geom_histogram(aes(x=IVW_meta_est,y=..density..))+
  geom_line(aes(x=IVW_meta_est,y=IVWmeta_N),color="red")+
  theme_Publication()+
  xlab("IVW meta Estimate")+
  ylim(limits=c(0,4))

p4 = ggplot(data1)+
  geom_histogram(aes(x=IVWs_est,y=..density..))+
  geom_line(aes(x=IVWs_est,y=IVWs_N),color="red")+
  theme_Publication()+
  xlab("IVW summary")+
  ylim(limits=c(0,4))


grid.arrange(p1,p2,p4,nrow=1)
ske_vec <- rep(0,4)
ske_vec_n <- rep(0,4)
kur_vec <- rep(0,4)
kur_vec_n <- rep(0,4)
for(k in 1:4){
  ske_vec[k] <- skewness(data1[,k])
  kur_vec[k] <- kurtosis(data1[,k])
  temp_n <- rnorm(10000,
                  mean = mean(data1[,k]),
                  sd = sd(data1[,k]))
  ske_vec_n[k] <- skewness(temp_n)
  kur_vec_n[k] <- kurtosis(temp_n)
}
result <- cbind(ske_vec,kur_vec,
                ske_vec_n,kur_vec_n)

row.names(result) <- c("Two Stage",
                       "IVW","IVW meta",
                       "IVW summary")

