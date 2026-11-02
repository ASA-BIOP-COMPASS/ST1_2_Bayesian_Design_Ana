# BayesianDesignAna: integer allocation and resource-saving search.
# Run after library(BayesianDesignAna).
ratio <- c(Placebo=1, Low=1, Mid=1, High=1)
rates <- c(Placebo=.30, Low=.40, Mid=.50, High=.60)

# 1:1:1:1 permits 164 or 168, not 165, 166, or 167:
allocation_unit(ratio)
valid_n <- candidate_totals(ratio, 120, 220)
allocate_arms(164, ratio)
allocate_arms(168, ratio)

# Fast preliminary search, then more simulations near the crossing:
ss <- design_sample_size(
  total_n=valid_n, ratios=ratio, response_rates=rates,
  borrowing_prop=.50, target_power=.80,
  search='two_stage', initial_step=5, refine_span=5,
  n_sim_initial=500, n_sim=10000, n_post=2000,
  multiplicity='hochberg', overall='any', seed=1234,
  table_max_cols=12, plot=TRUE)

ss$initial_results
ss$final_results
ss$minimum_candidate
plot(ss)   # No simulations repeated
