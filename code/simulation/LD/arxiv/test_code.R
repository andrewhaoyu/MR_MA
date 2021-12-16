Gamma_select = Gamma[idx.cau_m]
var_Gamma_select = var_Gamma[idx.cau_m]
alpha_select = alpha[idx.cau_m]
var_alpha_select = var_alpha[idx.cau_m]
0.07759045



Gamma = Gamma_select
se_Gamma = se_Gamma_select
alpha = alpha_select
se_alpha = se_alpha_select
ldscore = ld_score_select
MAF = MAF_select

dim(R)

alpha = alpha*sqrt(2*MAF*(1-MAF))
se_alpha=  se_alpha*sqrt(2*MAF*(1-MAF))
Gamma = Gamma*sqrt(2*MAF*(1-MAF))
se_Gamma = se_Gamma*sqrt(2*MAF*(1-MAF))
beta_est = as.numeric(crossprod(Gamma,alpha)/crossprod(alpha))
R_temp = R
R = diag(24)

diff_var = (se_Gamma^2+beta_est^2*se_alpha^2)
chi_est = (Gamma-alpha*beta_est)^2/diff_var
scale_ldscore = ldscore/diff_var
tau_est = coefficients(lm(chi_est~scale_ldscore-1))
tau_est = ifelse(tau_est>0,tau_est,0)

W = GetWtauMat(Gamma,se_Gamma,
               alpha,se_alpha,
               ldscore,tau_est,beta_est,R)
awa = quadform(x= as.matrix(alpha),M = W)
beta_est = awa^-1*
  crossprod(t(crossprod(alpha,W)),Gamma)

ObjFun <- function(Gamma,alpha,beta_est){
  W = GetWtauMat(Gamma,se_Gamma,
                 alpha,se_alpha,
                 ldscore,tau_est,beta_est,R)
  return(quadform(x= as.matrix(Gamma-beta_est*alpha),M = W))
}

ObjFun(Gamma,alpha,as.numeric(beta_est))
ObjFun(Gamma,alpha,as.numeric(0.05889707))

profile.loglike <- function(beta) {
  -(1/2) * sum((b_out - b_exp * beta)^2/(se_out^2 + se_exp^2 *
                                           beta^2))
}

b_out = Gamma_select
se_out = se_Gamma_select
b_exp = alpha_select
se_exp = se_alpha_select
bound <- quantile(abs(b_out/b_exp), 0.95, na.rm = TRUE) *
  2
beta.hat <- optimize(profile.loglike, bound * c(-1, 1), maximum = TRUE,
                     tol = .Machine$double.eps^0.5)$maximum
while (abs(beta.hat) > 0.95 * bound) {
  bound <- bound * 2
  beta.hat <- optimize(profile.loglike, bound * c(-1, 1),
                       maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
}


Gamma_select,se_Gamma_select,
                         alpha_select,se_alpha_select,
                         ld_score_select,R,MAF_select