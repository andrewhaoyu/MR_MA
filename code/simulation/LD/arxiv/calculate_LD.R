library(bigsnpr)
library(data.table)
library(dplyr)
cur.dir <- "/data/zhangh24/MR_MA/result/LD/"
j = 22
#snp_readBed(paste0(cur.dir,"chr",j,".hm3.bed"))
obj.bigSNP <- snp_attach(paste0(cur.dir,"chr",j,".hm3.rds"))
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
pos_table = fread(paste0(cur.dir,"../../data/ld_block"))
pos_chr = pos_table %>% 
  filter(chr==paste0("chr",j))
  
# parse LD blocks
# add the beginning and end to avoid dropping snps
if (pos_chr$stop[1] > 1) {
  pos_chr$start[1] = 1
}
pos_tail = data.frame(chr=paste0('chr',j), start=pos_chr$stop[length(pos_chr$stop)], stop=Inf)
pos_chr = rbind(pos_chr, pos_tail)


# cor(G[,idx],G[,jdx])
# sum(G[,idx])/(2*nrow(G))
# sum(G[,jdx])/(2*nrow(G))
NCORES = 2
n.sub = nrow(G)
#POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0(cur.dir), ncores = 2)
set.seed(666)
#size present window size for compute correlations
# corr0 <- snp_cor(G, ind.row = sample(c(1:n.sub),3000),
#                  ncores = NCORES,
#                  infos.pos = POS,
#                   size = 1000)
ix_list = list()
cor_list = list()


for (i in 1:nrow(pos_chr)) {
  cat(paste0(i, '..'))
  if (i == 1) {
    ix = which(CHR == j & POS >= pos_chr$start[i] & POS <= pos_chr$stop[i] )
  } else {
    ix = which(CHR == j & POS > pos_chr$start[i] & POS <= pos_chr$stop[i])
  }
  
  if (length(ix) == 0) {
    cor_list[[i]] = NULL
    ix_list[[i]] = NULL
    next
  } else {
    size = round((pos_chr$stop[i] - pos_chr$start[i]) / 1000)
    
    corr = snp_cor(G, ind.row = c(40001:80000),
                   ind.col=ix,
                  ncores = NCORES)
      # ev = eigen(corr)
      # cut = 286
      # if(nrow(corr)>cut){
      #   ev_mat = ev$vectors[,1:cut]
      #   eigen_value = ev$values[1:cut]
      #   # ev$vectors%*%diag(ev$values)%*%t(ev$vectors)
      #   corr_low = crossprod(t(ev_mat)*eigen_value,t(ev_mat))
      # }
      # 
      #snp_cor(G_tuning, ind.col=ix, size=size, alpha=1, thr_r2=0)
    
    saveRDS(corr, paste0("/data/zhangh24/MR_MA/data/eur_ld_block/chr", j, '_ld_', i, ".rds"))
    cor_list[[i]] = corr
    ix_list[[i]] = ix
  }
}
cat('\n')


library(Matrix)
M1 = matrix(rnorm(9),3,3)
M2 = matrix(rnorm(16),4,4)
M3 = matrix(rnorm(16),4,4)
corr0 = bdiag(cor_list)
save(corr0,file = paste0(cur.dir,"chr_",j,"_LDmat.rdata"))
# 
# library(rsvd)
# library(denoiseR)
# 
# 
# 
# 
# A = matrix(rnorm(9),3,3)
# B = A%*%t(A)
# 
# ev = eigen(B)
# ev$vectors%*%diag(ev$values)%*%t(ev$vectors)
# summary(ev2)
# ev
# ev2 = rrpca(B)
# adashrink(B)
# ev2$rotation%*%diag(ev2$eigvals)%*%ev2$rotation
# corr0 <- snp_cor(G, ind.row = c(40001:80000),
#                  ncores = NCORES,
#                  infos.pos = POS,
#                  size = 1000)
# 
# save(corr0,file = paste0(cur.dir,"chr_",j,"_LDmat.rdata"))
# G.temp = G[,1:1000]

load(paste0(cur.dir,"chr_",j,"_LDmat.rdata"))
save(G.temp,file = paste0(cur.dir,"example_G.rdata"))


res <- system(paste0("/data/zhangh24/software/plink2 --threads 2 --bfile /data/zhangh24/MR_MA/result/LD/chr",j,".sub.hm3 --out ",cur.dir,"chr_",j,"_LD --r bin"))
