function (b_exp, b_out, se_exp, se_out, diagnostics = FALSE)
{
  profile.loglike <- function(beta) {
    -(1/2) * sum((b_out - b_exp * beta)^2/(se_out^2 + se_exp^2 *
                                             beta^2))
  }
  bound <- quantile(abs(b_out/b_exp), 0.95, na.rm = TRUE) *
    2
  beta.hat <- optimize(profile.loglike, bound * c(-1, 1), maximum = TRUE,
                       tol = .Machine$double.eps^0.5)$maximum
  while (abs(beta.hat) > 0.95 * bound) {
    bound <- bound * 2
    beta.hat <- optimize(profile.loglike, bound * c(-1, 1),
                         maximum = TRUE, tol = .Machine$double.eps^0.5)$maximum
  }
  score.var <- sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 -
                                                         se_out^2) * se_exp^2 + se_exp^2 * se_out^2)/(se_out^2 +
                                                                                                        beta.hat^2 * se_exp^2)^2)
  I <- sum(((b_exp^2 - se_exp^2) * se_out^2 + (b_out^2 - se_out^2) *
              se_exp^2)/(se_out^2 + beta.hat^2 * se_exp^2)^2)
  dif <- b_out - beta.hat * b_exp
  dif.var <- se_out^2 + beta.hat^2 * se_exp^2
  chi.sq.test <- sum((dif/sqrt(dif.var))^2)
  if (diagnostics) {
    std.resid <- (b_out - b_exp * beta.hat)/sqrt((se_out^2 +
                                                    beta.hat^2 * se_exp^2))
    par(mfrow = c(1, 2))
    qqnorm(std.resid)
    abline(0, 1)
    beta.hat.loo <- rep(NA, length(b_out))
    beta.hat.loo <- rep(NA, length(b_out))
    if (length(b_out) > 100) {
      a <- quantile(abs(b_exp/se_exp), 1 - 100/length(b_out))
    }
    else {
      a <- 0
    }
    for (i in 1:length(b_out)) {
      if (abs(b_exp[i]/se_exp[i]) > a) {
        beta.hat.loo[i] <- mr.raps.simple(b_exp[-i],
                                          b_out[-i], se_exp[-i], se_out[-i])$beta.hat
      }
    }
    plot(abs(b_exp/se_exp), beta.hat.loo)
    abline(h = beta.hat)
    print(ad.test(std.resid))
    print(shapiro.test(std.resid))
  }
  out <- list(beta.hat = beta.hat, beta.se = sqrt(score.var/I^2),
              beta.p.value = min(1, 2 * (1 - pnorm(abs(beta.hat)/sqrt(score.var/I^2)))),
              naive.se = sqrt(1/I), chi.sq.test = chi.sq.test)
  if (diagnostics) {
    out$std.resid <- std.resid
    out$beta.hat.loo <- beta.hat.loo
  }
  out
}