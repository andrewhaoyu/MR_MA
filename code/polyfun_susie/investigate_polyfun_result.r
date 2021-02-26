#download data from https://alkesgroup.broadinstitute.org/polyfun_results
wget https://storage.googleapis.com/broad-alkesgroup-public/polyfun_results/body_BMIz.txt.gz
setwd("/data/zhangh24/MR_MA/")
data <- read.table(gzfile("./data/body_BMIz.txt.gz"),header=T)
dim(data)
idx <- which(data$pip>=0.95)
range(data$P_BOLT_LMM)
