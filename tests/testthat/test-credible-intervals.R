test_that('HPD intervals are shortest contiguous posterior-draw intervals', {
  x <- c(0.01, 0.03, 0.04, 0.06, 0.11, 0.15, 0.30, 0.55, 0.65, 0.90)
  h <- BayesianDesignAna:::.hpd_draw_interval(x, .80)
  expect_length(h, 2)
  expect_true(h[1] <= h[2])
  expect_true(sum(x >= h[1] & x <= h[2]) >= ceiling(.80 * length(x)))
  expect_equal(BayesianDesignAna:::.hpd_draw_interval(rep(.3, 10), .95), c(.3, .3))
})

test_that('all interval types share identical posterior simulation', {
  set.seed(1)
  dat <- data.frame(
    arm = rep(c('Placebo','Dose1','Dose2'), each=12),
    y = rep(c(0,1,1,0,0,1,0,1,0,1,0,1),3),
    stratum = rep(c('A','B','C'), 12)
  )
  cfg <- list(data=dat, covariates='stratum', n_iter=100, burnin=20,
              thin=2, n_chains=2, seed=789, multiplicity='none')
  both <- do.call(analyze_logistic, c(cfg, list(interval_method='both',
                                             credible_level=.90,
                                             sas_compatible_summary=TRUE)))
  hpd <- do.call(analyze_logistic, c(cfg, list(interval_method='HPD',
                                            credible_level=.90)))
  et <- do.call(analyze_logistic, c(cfg, list(interval_method='equal_tail',
                                           credible_level=.90)))
  expect_identical(both$posterior, hpd$posterior)
  expect_identical(both$posterior, et$posterior)
  expect_equal(both$summary$rd_mean, hpd$summary$rd_mean)
  expect_equal(both$summary$prob_superiority, et$summary$prob_superiority)
  expect_equal(both$summary$rd_lower, both$summary$equal_tail_lower)
  expect_equal(both$summary$rd_upper, both$summary$equal_tail_upper)
  expect_equal(hpd$summary$rd_lower, both$summary$hpd_lower)
  expect_equal(hpd$summary$rd_upper, both$summary$hpd_upper)
  expect_equal(et$summary$rd_lower, both$summary$equal_tail_lower)
  expect_equal(et$summary$rd_upper, both$summary$equal_tail_upper)
  expect_true(is.null(hpd$sas_summary))
  expect_identical(both$interval_method, 'both')
  expect_equal(both$credible_level, .90)
  expect_true(all(c('pP_ref_mean', 'pD1_ref_mean', 'pD2_ref_mean',
                    'RD_D1_mean', 'RD_D1_LCL', 'RD_D1_UCL',
                    'Pr_Sup_D1', 'Pr_Sup_D2') %in% colnames(both$sas_summary)))
  expect_equal(both$sas_summary$RD_D1_LCL, both$summary$hpd_lower[1])
  expect_equal(both$sas_summary$RD_D1_UCL, both$summary$hpd_upper[1])
  expect_equal(both$sas_summary$Pr_Sup_D1,
               mean(both$reference_prob_draws[, 'Dose1'] >
                    both$reference_prob_draws[, 'Placebo']))
})

test_that('interval arguments are validated', {
  dat <- data.frame(arm=rep(c('Placebo','Dose1'),each=4), y=rep(c(0,1),4))
  expect_error(analyze_logistic(dat,credible_level=1), 'credible_level')
  expect_error(analyze_logistic(dat,credible_level=0), 'credible_level')
  expect_error(analyze_logistic(dat,interval_method='not_a_method'), 'arg')
  expect_error(analyze_logistic(dat,sas_compatible_summary=NA),
               'sas_compatible_summary')
})
