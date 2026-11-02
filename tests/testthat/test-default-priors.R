test_that('default priors are reference/diffuse and no borrowing happens without opt in', {
  n <- c(Placebo=16, Low=18, High=18)
  rates <- c(Placebo=.3, Low=.4, High=.55)
  a <- design_power(n, rates, seed=77, n_sim=4, n_post=100)
  b <- design_power(n, rates, borrowed_n=0,
                    base_control_prior=c(.5,.5), treatment_prior=c(.5,.5),
                    seed=77, n_sim=4, n_post=100)
  expect_equal(a$overall_power, b$overall_power)
  expect_equal(a$scenario$borrow_ess, 0)
  expect_equal(unname(a$scenario$actual_control_prior), c(.5,.5))
  expect_equal(a$scenario$treatment_prior, c(.5,.5))
})

test_that('custom Beta priors flow through design_power and design_sample_size', {
  n <- c(Placebo=20, Dose1=20)
  rates <- c(Placebo=.3, Dose1=.5)
  fit <- design_power(n, rates, borrowed_n=10, historical_rate=.3,
                      base_control_prior=c(1,1), treatment_prior=c(2,2),
                      seed=15,n_sim=3,n_post=100)
  expect_equal(unname(fit$scenario$actual_control_prior), c(4,8))
  expect_equal(fit$scenario$treatment_prior, c(2,2))

  ss <- design_sample_size(total_n=c(40,60),ratios=c(Placebo=1,Dose1=1),
                          response_rates=rates, borrowed_n=0,
                          base_control_prior=c(1,1),treatment_prior=c(2,2),
                          seed=9,n_sim=3,n_post=100)
  expect_equal(ss$base_control_prior, c(1,1))
  expect_equal(ss$treatment_prior, c(2,2))
})

test_that('invalid prior shapes fail with informative errors', {
  n <- c(Placebo=16, Dose1=16)
  rates <- c(Placebo=.3,Dose1=.5)
  expect_error(design_power(n, rates, base_control_prior=c(0,1),n_sim=2,n_post=100),
               'base_control_prior')
  expect_error(design_power(n, rates, treatment_prior=c(1,-2),n_sim=2,n_post=100),
               'treatment_prior')
})

test_that('logistic default has no historical placebo prior', {
  set.seed(19)
  dat <- data.frame(arm=rep(c('Placebo','Dose1'),each=15),
                    severity=factor(rep(c('Mild','Severe'),15)),
                    y=rep(c(0,1,0,1,0),6))
  fit <- analyze_logistic(dat,covariates='severity',
                          seed=19,n_iter=120,burnin=40,thin=2,n_chains=2)
  expect_equal(fit$control_prior, c(.5,.5))
  expect_equal(fit$intercept_sd, 100)
  expect_equal(fit$covariate_sd, 100)

  custom <- analyze_logistic(dat,covariates='severity',
                             control_prior=c(5.52,14.48),
                             intercept_sd=2.5,covariate_sd=1.5,
                             seed=19,n_iter=120,burnin=40,thin=2,n_chains=2)
  expect_equal(custom$control_prior, c(5.52,14.48))
  expect_equal(custom$intercept_sd, 2.5)
  expect_equal(custom$covariate_sd, 1.5)
})

test_that('multiple scenarios propagate custom priors', {
  d <- matrix(c(.1,.2), ncol=2,dimnames=list(NULL,c('Low','High')))
  res <- plot_power(delta_matrix=d, total_n=c(60,90),
                    ratios=c(Placebo=1,Low=1,High=1),placebo_rate=.3,
                    plot=FALSE,base_control_prior=c(1,1),
                    treatment_prior=c(2,2),seed=20,n_sim=3,n_post=100)
  expect_equal(res$scenarios[[1]]$base_control_prior,c(1,1))
  expect_equal(res$scenarios[[1]]$treatment_prior,c(2,2))
})
