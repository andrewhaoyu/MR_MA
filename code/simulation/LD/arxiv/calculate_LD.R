library(bigsnpr)
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
j = 22
#snp_readBed(paste0(cur.dir,"chr",j,".hm3.bed"))
obj.bigSNP <- snp_attach(paste0(cur.dir,"chr",j,".hm3.rds"))
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
NCORES = 2
n.sub = nrow(G)
#POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0(cur.dir), ncores = 2)
set.seed(666)
#size present window size for compute correlations
corr0 <- snp_cor(G, ind.row = sample(c(1:n.sub),3000),
                 ncores = NCORES,
                 infos.pos = POS,
                  size = 1000)
save(corr0,file = paste0(cur.dir,"chr_",j,"_LDmat.rdata"))
G.temp = G[,1:1000]

load(paste0(cur.dir,"chr_",j,"_LDmat.rdata"))
save(G.temp,file = paste0(cur.dir,"example_G.rdata"))


res <- system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /data/zhangh24/MR_MA/result/LD/chr",j,".sub.hm3 --out ",cur.dir,"chr_",j,"_LD --r bin"))
