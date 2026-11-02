test_that("overall any is default and all is optional", {
  n <- c(Placebo=20,Low=20,High=20)
  r <- c(Placebo=.3,Low=.4,High=.6)
  a <- design_power(n,r,borrowed_n=5,n_sim=5,n_post=150)
  b <- design_power(n,r,borrowed_n=5,overall="all",n_sim=5,n_post=150)
  expect_equal(a$overall_definition,"any")
  expect_equal(b$overall_definition,"all")
  expect_gte(a$overall_power,b$overall_power)
})
test_that("single and multi scenario figures can be created", {
  ratios <- c(Placebo=1,Low=2,Mid=2,High=1)
  rates <- c(Placebo=.3,Low=.4,Mid=.5,High=.6)
  ss <- design_sample_size(c(60,120,180),ratios,rates,
      borrowing_prop=.3, n_sim=3,n_post=120)
  grDevices::pdf(file=tempfile(fileext=".pdf"),width=15,height=9)
  on.exit(grDevices::dev.off(),add=TRUE)
  expect_silent(plot(ss))
  mat <- rbind(c(.05,.10,.20),c(.1,.2,.3))
  colnames(mat) <- names(ratios)[-1]
  expect_silent(plot_power(mat,c(60,120,180),ratios,.3,
    borrowing_prop=.3,n_sim=3,n_post=120))
})
