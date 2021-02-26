#
args = commandArgs(trailingOnly = T)
i1 = as.numeric(args[[1]])
# library(devtools)
# install_github("andrewhaoyu/bc2")

setwd("/data/zhangh24/UKB_LD/")
library(bc2)
num <- 1000

ld_filename <- read.table("ukb_ld_file_names.txt",header=F)
n <- nrow(ld_filename)
start.end <- startend(n,num,i1)
start = start.end[1]
end = start.end[2]
for(l in start:end){
  system(paste0("cd /data/zhangh24/UKB_LD/ | wget https://storage.googleapis.com/broad-alkesgroup-public/UKBB_LD/",ld_filename[l,1]))  
}

