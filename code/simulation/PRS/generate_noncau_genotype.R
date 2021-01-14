arg <- commandArgs(trailingOnly=T)
#i1 = 2,..30 generate 29 different genotype datasets for non causal snps
i1 <- as.numeric(arg[[1]])

#generate independent causal SNPs
set.seed(i1)
n.sub <- 110000
n.snp <- 5000

MAF <- 0.25

genotype <- matrix(rbinom(n.sub*n.snp,2,0.25),n.sub,n.snp)
#genotype of M
genotype_s <- (genotype-2*MAF)/(sqrt(2*MAF*(1-MAF)))
save(genotype_s,file = paste0("/data/zhangh24/MR_MA/result/simulation/prs/noncau_genotype_M_",i1,".rdata"))
n.sub <- 100000
n.snp <- 5000
genotype2 <- matrix(rbinom(n.sub*n.snp,2,0.25),n.sub,n.snp)
#genotype of Y
genotype_s2 <- (genotype2-2*MAF)/(sqrt(2*MAF*(1-MAF)))
save(genotype_s2,file = paste0("/data/zhangh24/MR_MA/result/simulation/prs/noncau_genotype_Y_",i1,".rdata"))


