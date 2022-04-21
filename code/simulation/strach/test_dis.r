#test the distribution of spike and slab
p = 0.3
n = 100000
eta = rbinom(n,1,p)
idx <- which(eta!=0)
beta = rep(0,n)
beta[idx] = rnorm(length(idx),0,1)
hist(beta)
quantile(beta,probs = c(0.03,0.05,0.10,0.15,0.5,0.85,0.90,0.95,0.97))


beta2 = rbinom(n,1,p)*rnorm(n,0,1)
quantile(beta2,probs = c(0.03,0.05,0.10,0.15,0.5,0.85,0.90,0.95,0.97))
