#generate independent causal SNPs
n.sub <- 110000
n.snp <- 5000

MAF <- 0.25

genotype <- matrix(rbinom(n.sub*n.snp,2,0.25),n.sub,n.snp)
#genotype of M
genotype_s <- (genotype-2*MAF)/(sqrt(2*MAF*(1-MAF)))
save(genotype_s,file = "/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_M.rdata")
n.sub <- 100000
n.snp <- 5000
genotype2 <- matrix(rbinom(n.sub*n.snp,2,0.25),n.sub,n.snp)
#genotype of Y
genotype_s2 <- (genotype2-2*MAF)/(sqrt(2*MAF*(1-MAF)))
save(genotype_s2,file = "/data/zhangh24/MR_MA/result/simulation/prs/cau_genotype_Y.rdata")


