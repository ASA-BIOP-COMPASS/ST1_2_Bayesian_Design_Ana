test_that('allocation sums with unequal ratios', {
  x <- allocate_arms(102, c(Control = 1, A = 2, B = 3))
  expect_equal(sum(x), 102)
  expect_length(x, 3)
})
test_that('borrowing caps respect constraints', {
  x <- borrowing_limits(20, 100, max_fraction = .5, max_absolute = 12)
  expect_equal(x$allowed_ess, 12)
  expect_equal(borrowing_limits(0, 50)$allowed_ess, 0)
})
test_that('two and three treatment arms can be analyzed', {
  x <- analyze_binary(c(Control = 4, A = 8, B = 10),
                      c(Control = 12, A = 15, B = 15), n_post = 200, seed = 1)
  expect_equal(nrow(x), 2)
  expect_true(all(x$prob_superiority >= 0 & x$prob_superiority <= 1))
})
test_that('simulation and sample size search produce outputs', {
  r <- c(Control = .3, A = .5, B = .55)
  a <- simulate_power(c(Control = 10, A = 12, B = 12), r, n_sim = 5, n_post = 100)
  expect_equal(nrow(a$per_arm), 2)
  b <- find_sample_size(c(21, 30), c(Control = 1, A = 1, B = 1), r,
                        n_sim = 5, n_post = 100)
  expect_equal(nrow(b$results), 2)
})
test_that('logistic analysis supports arbitrary arm names', {
  set.seed(2)
  d <- data.frame(arm = rep(c('Placebo', 'A', 'B'), each = 12),
                  strat1 = rep(c(0, 1), 18), y = rbinom(36, 1, .4))
  fit <- analyze_logistic(d, covariates = 'strat1', n_iter = 150,
                          burnin = 50, thin = 2, n_chains = 2)
  expect_equal(nrow(fit$summary), 2)
  expect_true(is.finite(fit$summary$rd_mean[1]))
})
