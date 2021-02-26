setwd("/data/zhangh24/MR_MA")
library(data.table)
library(dplyr)
library(tidyr)


#set up the python environment for polyfun using install_software.R
#polyfun needs CHR, BP, A1, A2
bmi = as.data.frame(fread("./data/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt",header=T))
n = nrow(bmi)
idx <- which(duplicated(bmi$SNP))
bmi = bmi %>% 
  mutate(A1 = Tested_Allele,
         A2 = Other_Allele,
         BP = POS,
         Z = BETA/SE,
         N = N) %>% 
  select(SNP,CHR,BP,A1,A2,BETA,Z,N,P) %>% 
  distinct(SNP,.keep_all=T)

idx <- which(duplicated(bmi$SNP))

write.table(bmi,file = paste0("./result/real_data_analysis/bmi/bmi_sumdata.txt"),col.names = T,row.names = F,quote=F)
cur.dir <- "/data/zhangh24/MR_MA"
bmi.filter <- bmi %>% 
  filter(CHR==1) %>% 
  filter(BP>=1708801-3E-06&
           BP<= 753541+3E+06)
#use polyfun to compute prior causal probability
#paste0("python extract_snpvar.py --sumstats ",cur.dir,"/result/real_data_analysis/bmi/bmi_sumdata.txt --out ",cur.dir,"/result/real_data_analysis/bmi/bmi_polyfunout --allow-missing")
#paste0("python polyfun.py --compute-h2-L2 --no-partitions --output-prefix output/testrun --sumstats cur.dir,"/result/real_data_analysis/bmi/bmi_sumdata.txt --out ",cur.dir,"/result/real_data_analysis/bmi/bmi_polyfunout --ref-ld-chr example_data/annotations. --w-ld-chr example_data/weights.")

python finemapper.py 
--geno /data/zhangh24/KG.plink/EUR 
--sumstats /data/zhangh24/MR_MA/result/real_data_analysis/bmi/bmi_polyfunout
--chr 1 \
--start 753541
--end 3753541
--method susie \
--max-num-causal 5 \
--cache-dir LD_cache \
--out output/test






output <- as.data.frame(fread("/gpfs/gsfs11/users/zhangh24/MR_MA/software/polyfun/output/test",header=T))

idx <- which(bmi.filter=="rs6603803")
bmi.filter[idx,]

library(susieR)
set.seed(1)
n    <- 1000
p    <- 1000
beta <- rep(0,p)
beta[c(1,2,300,400)] <- 1
X   <- matrix(rnorm(n*p),nrow=n,ncol=p)
y   <- X %*% beta + rnorm(n)
res <- susie(X,y,L=10)
plot(coef(res),pch = 20)

N = 600
P = 1000

Sigma = matrix(0,P,P)
diag(Sigma) = 1
Sigma[upper.tri(Sigma)] = 0.3
Sigma[lower.tri(Sigma)] = 0.3

library(mvtnorm)
X = rmvnorm(N,mean = rep(0,P),sigma = Sigma)
beta = rep(0,P)
beta[c(1,10,300)] = 1
Y = X%*%beta + rnorm(N)
res <- susie(X,Y,L=10)
coef(res)[301]
sumstats = univariate_regression(X,Y)
z_scores <- sumstats$betahat / sumstats$sebetahat
susie_plot(z_scores, y = "z", b=beta)
fitted <- susie(X, Y,
                L = 10,
                estimate_residual_variance = TRUE, 
                estimate_prior_variance = FALSE,
                scaled_prior_variance = 0.1,
                verbose = TRUE)
print(fitted$sets)
susie_plot(fitted, y="PIP", b=beta)
i  <- fitted$sets$cs[[3]]
z3 <- cbind(i,z_scores[i],fitted$pip[i])
colnames(z3) <- c('position', 'z-score', 'PIP')
z3[order(z3[,2], decreasing = TRUE),]

R = cor(X)
fitted_rss <- susie_rss(z_scores, R, L = 10,
                        estimate_residual_variance = TRUE, 
                        estimate_prior_variance = TRUE)

susie_plot(fitted_rss, y="PIP", b=beta)
