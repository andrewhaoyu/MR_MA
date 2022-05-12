function (data, diagnostics = TRUE, over.dispersion = TRUE, loss.function = "huber",
          ...)
{
  prior.param <- fit.mixture.model(data$beta.exposure/data$se.exposure)
  out <- mr.raps.shrinkage(data$beta.exposure, data$beta.outcome,
                           data$se.exposure, data$se.outcome, over.dispersion, loss.function,
                           prior.param = prior.param, diagnostics = diagnostics,
                           ...)
  if (diagnostics) {
    cat(paste0("Estimated causal effect: ", signif(out$beta.hat,
                                                   3), ", standard error: ", signif(out$beta.se, 3),
               ", p-value: ", signif(pnorm(-abs(out$beta.hat/out$beta.se)) *
                                       2, 3), ".\n"))
    cat(paste0("Estimated overdispersion variance: ", signif(out$tau2.hat,
                                                             3), ", standard error: ", signif(out$tau2.se, 3),
               ", p-value: ", signif(pnorm(-abs(out$tau2.hat/out$tau2.se)) *
                                       2, 3), ".\n"))
    cat(paste0("ANOVA test: are the weights and residuals independent? \n"))
    weights <- out$gamma.hat.z
    std.resids <- out$t
    df <- max(round(length(weights)/20), 3)
    lm.test <- lm(std.resids ~ bs(weights, knots = quantile(weights,
                                                            1:df/(df + 1))) - 1)
    print(anova(lm.test))
  }
  out
}