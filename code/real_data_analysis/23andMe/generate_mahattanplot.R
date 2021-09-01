TransArraytoVec <- function(array,vec){
  n = length(vec)
  command = paste0("expand.grid(")
  for(k in 1:(n-1)){
    temp = paste0("c(1:",vec[k],"),")
    command = paste0(command,temp)
  }
  temp = paste0("c(1:",vec[n],"))")
  command = paste0(command,temp)
  result =  eval(parse(text = command))
  return(result[array,])
}
setwd("/n/holystore01/LABS/xlin/Lab/hzhang/MR_MA/")
library(data.table)
library(R.utils)
library(tidyverse)
library(qqman)
library(dplyr)

args = commandArgs(trailingOnly = T)
array = as.numeric(args[[1]])
vecargs = TransArraytoVec(array,c(4,5))
i1 = as.numeric(vecargs[1])
l = as.numeric(vecargs[2])
eth = c("european","african_american",
        "latino","trans_ethnic")
ethname = c("European",
            "African American",
            "Latino",
            "Trans ethnic")
trait = c("positive_vs_negative",
           "positive_hospitalized_dx_negative_controls",
           "positive_pneumonia_dx_negative_controls",
           "positive_respiratory_broad_dx_negative_controls",
           "positive_respiratory_support_dx_negative_controls")
traitname = c("Test-positive versus test-negative",
               "Hospitalized",
               "Pneumonia",
               "Respriatory support",
               "Severe respiratory sympotoms"
)


load(paste0("./data/cleaned/",eth[i1],"/",trait[l],".rdata"))
data = data %>% 
  rename(BP = position,
         P = pvalue) %>% 
  select(rsid,CHR,BP,P,MAF)
dat = data
sample_size <- as.data.frame(fread("./data/cleaned/23andme_sample_size.csv",header=T))
library(readr)
x = dat$P
z = qnorm(x / 2)
lambda = round(median(z^2) / qchisq(0.5,1), 3)

idx <- which(sample_size$eth==ethname[i1]&
               sample_size$Disease==traitname[l])
N.effect  <- as.integer(sample_size[idx,"N_effect"])
#rescale lambda to 1000 subjects

lambda_1000 = round(1+500*(lambda-1)/N.effect ,3)

convert.qval.pval = function(qvalues) {
  # you need to know the estimate of pi0 used to create the q-value
  # that's the maximum q-value (or very, very close to it)
  pi0 = max(qvalues)
  # compute m0, the estimated number of true nulls
  m0 = length(qvalues) * pi0
  # then you multiply each q-value by the proportion of true nulls
  # expected to be under it (the inverse of how you get there from
  # the p-value):
  return(qvalues * rank(qvalues) / m0)
}
p.pwas <- 5E-08
#q.pwas <- convert.qval.pval(c(tmp, 0.05))[(nrow(dat)+1)]

nCHR <- length(unique(dat$CHR))
dat$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(dat$CHR)){
  nbp[i] <- max(dat[dat$CHR == i,]$BP)
  dat$BPcum[dat$CHR == i] <- dat$BP[dat$CHR == i] + s
  s <- s + nbp[i]
}
library(dplyr)
axis.set <- dat %>% 
  group_by(CHR) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2)
ylim <- abs(floor(log10(min(dat$P)))) + 2 
sig1 <- p.pwas
#sig2 <- q.pwas

#sigline <- data.frame(sig=c(-log10(sig1),-log10(sig2)),val=c(paste0("P=",signif(sig1,2)),"FDR=0.05"))
sigline <- data.frame(sig=c(-log10(sig1)),val=c(paste0("P=",signif(sig1,2))))
library(ggplot2)
manhplot <- ggplot(dat, aes(x = BPcum, y = -log10(P), 
                            color = as.factor(CHR), size = -log10(P))) +
  geom_point(alpha = 0.8, size=0.8) + 
scale_x_continuous(label = axis.set$CHR, breaks = axis.set$center) +
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
  scale_color_manual(values = rep(c("#08306b", "#4292c6"), nCHR)) +
  scale_size_continuous(range = c(0.5,3)) +
  geom_hline(data = sigline, aes(yintercept = sig), color= "red", linetype="dashed") +
  guides(color = F) + 
  labs(x = NULL, 
       y = "-log10(p)", 
       linetype = "",
       title = paste0(traitname[l]," for ",ethname[i1]))+
  #subtitle = "A2: Critically ill COVID19+ vs. population controls;\nB1: Hospitalized COVID19+ vs non-hospitalized COVID19+;\nB2: Hospitalized COVID19+ vs. population controls;\nC2: Reported SARS-CoV-2 infection vs. population controls") + 
  theme(
    legend.position = "top",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5),
    plot.title = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 8)
  )
outpath = "/n/holystore01/LABS/xlin/Lab/hzhang/MR_MA/result/23andMe/qqman/"
ggsave(filename=paste0("man_",eth[i1],"_",trait[l],".png"),
       plot=manhplot, device="png",
       path=outpath,
       width=9, height=4, units="in", dpi=300)


png(filename = paste0(outpath,"/QQ_",eth[i1],"_",trait[l],".png"), width = 8, height = 8, units = "in",res=300)
qq(dat$P)
text(5,1,expression(paste(lambda[1000]," = ")),cex = 1.5)
text(5.7,1,paste(lambda_1000),cex = 1.5)
title(paste0(trait_name[l]," for ",eth[i1]))
dev.off()

