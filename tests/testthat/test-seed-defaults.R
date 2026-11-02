test_that("Hochberg default and seed can be controlled", {
  n <- c(Placebo=20,Low=20,High=20)
  rates <- c(Placebo=.3,Low=.4,High=.6)
  a <- design_power(n, rates, borrowed_n=5, seed=42, n_sim=6, n_post=150)
  b <- design_power(n, rates, borrowed_n=5, seed=42, n_sim=6, n_post=150)
  expect_equal(a$overall_power, b$overall_power)
  expect_equal(a$per_arm, b$per_arm)
  expect_equal(a$scenario$multiplicity, "hochberg")
  ratios <- c(Placebo=1, Low=1, High=1)
  ss <- design_sample_size(c(60,90), ratios, rates, seed=42,
                           borrowed_n=5, n_sim=4,n_post=100)
  expect_equal(ss$seed,42)
  expect_equal(ss$multiplicity,"hochberg")
})
