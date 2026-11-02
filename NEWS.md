# BayesianDesignAna 0.6.0

- Added `interval_method = 'both'`, `'HPD'`, or `'equal_tail'` and user-selected `credible_level` to `analyze_logistic()`.
- Added an opt-in one-row `sas_summary` with posterior placebo/dose means, risk differences, HPD limits and superiority probabilities.
- Kept `rd_lower` and `rd_upper` equal-tailed when `interval_method='both'` to preserve previous behavior.
- Renamed the repository references to `JinhaoZou/ST1_1_Bayesian_Design_Ana`; installed R package stays `BayesianDesignAna`.
- Updated README with standalone examples, installation instructions, references and statistical limitations.

# BayesianDesignAna 0.5.0

- Exact integer-ratio allocation: 1:1:1:1 totals must be multiples of four; off-ratio totals never receive rounded allocations.
- New `allocation_unit()` and `candidate_totals()` helpers generate permissible blocks and inclusive feasible ranges.
- `design_sample_size(search="two_stage")` now screens with `n_sim_initial`, then reruns only a small window with `n_sim` high-precision trials and selects only among final-stage results.
- `results` stores combined initial/final points, `initial_results` and `final_results` preserve auditability, and `stage` plus `n_sim` identify each estimate's precision.
- Reworked single/multiple scenario tables and crossing labels to use stored results and display manageable subsets of candidate totals.
- Replaced README with a complete API/options reference and a runnable short example for every exported function.
- Does not change Bayesian prior defaults or statistical decision rules.
- Research preview; R CMD check and statistical validation remain required.

# BayesianDesignAna 0.4.3

- Changed default `analyze_logistic()` placebo prior to Beta(0.5, 0.5), replacing the example-specific Beta(5.52, 14.48). Historical information is never incorporated automatically.
- Kept diffuse proper Normal(0,100^2) defaults for active-dose log-odds intercepts and stratification effects; users can set `intercept_sd` and `covariate_sd` explicitly.
- Added explicit `base_control_prior` and `treatment_prior` options to `simulate_power()`, `design_power()`, `design_sample_size()`, and `plot_power()`, defaulting to Beta(0.5,0.5).
- Records the actual effective Beta prior for the placebo in `simulate_power()$scenario$actual_control_prior`; sample-size objects record baseline prior choices.
- Documents reference priors and optional historical borrowing; expands tests.
- Preserved existing power, multiple-arm, borrowing selection, and plotting behavior.

Release validation and reproducibility checks are required before scientific use.

