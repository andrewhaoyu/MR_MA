library(bigsnpr)
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
j = 22
#snp_readBed(paste0(cur.dir,"chr",j,".hm3.bed"))
obj.bigSNP <- snp_attach(paste0(cur.dir,"chr",j,".hm3.rds"))
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
NCORES = 2
n.snp = nrow(G)
#POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0(cur.dir), ncores = 2)
set.seed(666)
corr0 <- snp_cor(G, ind.row = sample(c(1:n.snp),500),
                 ncores = NCORES,
                 infos.pos = POS,
                  size = 500)
save(corr0,file = paste0(cur.dir,"chr_",j,"_LDmat.rdata"))
