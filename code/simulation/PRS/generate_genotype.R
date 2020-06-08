arg <- commandArgs(trailingOnly=T)
i1 <- as.numeric(arg[[1]])

#generate independent causal SNPs
library(bigmemory)
n.sub <- 110000
n.snp <- 200

MAF <- 0.25
genotype <- matrix(rbinom(n.sub*n.snp,2,0.25),n.sub,n.snp)
genotype_s <- (genotype-2*MAF)/(2*sqrt(MAF*(1-MAF)))

n.sub <- 100000
n.snp <- 150000
genotype2 <- matrix(rbinom(n.sub*n.snp,2,0.25),n.sub,n.snp)
genotype_s2 <- (genotype2-2*MAF)/(2*sqrt(MAF*(1-MAF)))
