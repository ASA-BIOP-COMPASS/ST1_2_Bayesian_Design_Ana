test_that("logistic analysis supports multiple multilevel factors", {
  n <- c(Placebo=12, Low=12, Mid=12, High=12)
  rates <- c(Placebo=.30, Low=.40, Mid=.50, High=.60)
  dat <- simulate_logistic_data(n, rates,
    strata=list(region=c(North=.4,South=.35,West=.25),
                severity=c(Mild=.65,Moderate=.25,Severe=.10)),
    stratum_log_or=list(region=c(.15,-.20),severity=c(.30,.55)))
  fit <- analyze_logistic(dat, covariates=c("region","severity"),
                          n_iter=80,burnin=20,n_chains=2,seed=12)
  expect_equal(ncol(fit$posterior), 8L)
  expect_equal(length(unique(fit$mcmc$parameter)),8L)
  expect_true(all(c("beta_regionSouth", "beta_regionWest",
    "beta_severityModerate", "beta_severitySevere") %in%
    colnames(fit$posterior)))
})
