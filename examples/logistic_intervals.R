# Credible intervals and SAS-style summary from the same posterior draws.
library(BayesianDesignAna)

n <- c(Placebo=40, Dose1=40, Dose2=40)
rates <- c(Placebo=.30, Dose1=.45, Dose2=.55)
dat <- simulate_logistic_data(n, rates,
  strata=list(severity=c(Mild=.65, Severe=.35)),
  stratum_log_or=list(severity=.30), seed=2026)

fit <- analyze_logistic(dat, covariates='severity', center='overall',
  control_prior=c(.5,.5), interval_method='both', credible_level=.95,
  sas_compatible_summary=TRUE, multiplicity='none',
  n_iter=5000,burnin=1000,n_chains=2,seed=2026)
fit$summary
fit$sas_summary
fit$mcmc
