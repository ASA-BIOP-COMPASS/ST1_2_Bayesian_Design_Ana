test_that('two borrowing modes and ESS reference work', {
 n <- c(Placebo=10,Dose=10); rates <- c(Placebo=.3,Dose=.5)
 expect_equal(design_power(n,rates,borrowed_n=5,n_sim=2,n_post=100)$borrowed_ess,5)
 expect_equal(borrowing_ess(30,.4),20)
 expect_error(design_power(n,rates,borrowed_n=2,borrowing_prop=.2,n_sim=2,n_post=100))
 ss <- design_sample_size(c(20,40),c(Placebo=1,Dose=1),rates,
     borrowed_n=8,historical_ess=5,n_sim=2,n_post=100)
 expect_true(all(ss$results$exceeds_historical_ess))
 expect_s3_class(ss,'bda_sample_size')
})
test_that('dose matrix curves and multi-level stratification', {
 ratio <- c(Placebo=1,Low=2,Mid=2,High=1)
 dm <- rbind(c(.1,.2,.3),c(.05,.15,.25))
 colnames(dm) <- names(ratio)[-1]
 out <- plot_power(dm,c(60,120),ratio,.3,borrowed_n=5,
   n_sim=2,n_post=100,plot=FALSE)
 expect_equal(nrow(out$results),4)
 dat <- simulate_logistic_data(c(Placebo=20,Low=20),
   c(Placebo=.3,Low=.5), strata=list(region=c(A=.2,B=.3,C=.5)))
 expect_equal(nlevels(dat$region),3)
 fit <- analyze_logistic(dat,covariates='region',n_iter=40,burnin=10,n_chains=2)
 expect_equal(nrow(fit$summary),1)
})
