test_that('exact 1:1:1:1 and 1:2:2:1 allocation blocks are enforced', {
  r <- c(Placebo=1,Low=1,Mid=1,High=1)
  expect_equal(allocation_unit(r),4)
  expect_equal(candidate_totals(r,160,172),c(160L,164L,168L,172L))
  expect_equal(unname(allocate_arms(164,r)),rep(41L,4))
  expect_equal(unname(allocate_arms(168,r)),rep(42L,4))
  expect_error(allocate_arms(165,r),'multiple of 4')
  expect_error(allocate_arms(166,r),'multiple of 4')
  expect_equal(allocation_unit(c(Placebo=2,Low=2,Mid=2,High=2)),4)
  expect_equal(allocation_unit(c(Placebo=1,Low=2,Mid=2,High=1)),6)
  expect_equal(unname(allocate_arms(126,c(Placebo=1,Low=2,Mid=2,High=1))),
               c(21L,42L,42L,21L))
  expect_equal(candidate_totals(c(Placebo=1,Low=2,Mid=2,High=1),120,144),
               c(120L,126L,132L,138L,144L))
})

test_that('candidate_totals honors boundaries and block stride', {
  r <- c(Placebo=1,Dose=1)
  expect_equal(candidate_totals(r,10,20,step_blocks=2),c(10L,14L,18L))
  expect_equal(candidate_totals(r,11,11),integer())
  expect_error(candidate_totals(r,12,20,step_blocks=0),'positive integer')
  expect_error(allocation_unit(c(Placebo=1.5,Dose=2)),'integer weights')
})

test_that('grid drops impossible totals instead of rounding arms', {
  ratio <- c(Placebo=1,Low=1,Mid=1,High=1)
  rates <- c(Placebo=.3,Low=.4,Mid=.5,High=.6)
  expect_warning(
    out <- design_sample_size(160:169,ratio,rates,
      n_sim=2,n_post=100,seed=8), 'Skipped')
  expect_equal(out$results$total_n,c(160L,164L,168L))
  expect_equal(out$results$Placebo,c(40L,41L,42L))
  expect_equal(out$results$High,c(40L,41L,42L))
  expect_equal(out$search,'grid')
})

test_that('two-stage search evaluates high precision only near crossing', {
  ratio <- c(Placebo=1,Low=1,Mid=1,High=1)
  rates <- c(Placebo=.3,Low=.4,Mid=.5,High=.6)
  grid <- candidate_totals(ratio,120,180)
  ss <- design_sample_size(grid,ratio,rates,
    search='two_stage',initial_step=4,refine_span=4,
    n_sim_initial=4,n_sim=20,n_post=100,seed=71)
  expect_s3_class(ss,'bda_sample_size')
  expect_true(all(ss$results$total_n %% 4 == 0))
  expect_equal(ss$allocation_block,4)
  expect_true(all(ss$initial_results$n_sim==4L))
  expect_true(all(ss$final_results$n_sim==20L))
  expect_true(all(ss$final_results$stage=='final'))
  expect_equal(anyDuplicated(ss$results$total_n),0L)
  expect_true(nrow(ss$initial_results) < length(grid))
  expect_true(nrow(ss$final_results) < length(grid))
  expect_equal(ss$results$overall_power[match(ss$final_results$total_n,
                                                  ss$results$total_n)],
               ss$final_results$overall_power)
  if (!is.null(ss$minimum_candidate)) {
    expect_true(ss$minimum_candidate$total_n %in% ss$final_results$total_n)
    expect_equal(ss$minimum_candidate$n_sim,20L)
  }
})

test_that('plotting never changes the cached simulation results', {
  ratio <- c(Placebo=1,Dose=1)
  rates <- c(Placebo=.3,Dose=.6)
  ss <- design_sample_size(candidate_totals(ratio,40,60),ratio,rates,
    search='two_stage',n_sim_initial=3,n_sim=8,n_post=100,
    initial_step=3,refine_span=3,seed=30,table_max_cols=4)
  before <- ss$results
  grDevices::pdf(file=tempfile(fileext='.pdf'),width=12,height=7)
  on.exit(grDevices::dev.off(),add=TRUE)
  expect_silent(graphics::plot(ss))
  expect_equal(before,ss$results)
})

test_that('plot_power supports two-stage different scenario grids', {
  ratio <- c(Placebo=1,Low=1,High=1)
  m <- rbind(c(.08,.12),c(.15,.25))
  colnames(m) <- c('Low','High')
  grDevices::pdf(file=tempfile(fileext='.pdf'),width=14,height=8)
  on.exit(grDevices::dev.off(),add=TRUE)
  expect_silent(out <- plot_power(m,candidate_totals(ratio,60,90),ratio,.30,
    search='two_stage',initial_step=3,refine_span=3,
    n_sim_initial=3,n_sim=8,n_post=100,seed=9,table_max_cols=5))
  expect_length(out$scenarios,2)
  expect_true(all(out$results$total_n %% 3 == 0))
})
