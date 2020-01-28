setwd("/Users/zhangh24/GoogleDrive/MR_MA")
#load CAD data
library(data.table)

#load LDL data
LDL <- as.data.frame(fread("./data/jointGwasMc_LDL.txt",header=T))
colnames(LDL)[9] <- "P"
idx <- which(as.numeric(LDL$P)<=5E-08)
length(idx)
