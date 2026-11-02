# BayesianDesignAna

**Bayesian sample-size planning and analysis for multi-arm trials with binary endpoints**  
R package: `BayesianDesignAna` · Version: **0.6.0** · GitHub repository: **`ST1_2_Bayesian_Design_Ana`**

`BayesianDesignAna` helps researchers plan trials with a shared placebo/control arm, compare multiple doses, and assess whether historical control information could be borrowed. It includes simple Beta-binomial simulations and a separate patient-level Bayesian logistic regression model with covariate adjustment. Outputs are intended for research and method development, not unvalidated confirmatory decisions.

## What can I do with this package?

| Goal | Main function |
| --- | --- |
| Generate sample sizes that follow an exact randomization ratio | `candidate_totals()` and `allocate_arms()` |
| Estimate power for a specified sample size | `design_power()` |
| Find a sample size that reaches a target power | `design_sample_size()` |
| Compare power curves across different dose effects | `plot_power()` |
| Specify and limit historical-control borrowing | `borrowing_ess()`, `borrowing_prior()`, `borrowing_limits()` |
| Estimate historical MAP prior effective sample size | `historical_map_ess()` (uses **RBesT**) |
| Analyze observed responder counts | `analyze_binary()` |
| Analyze patient-level data with stratification factors | `analyze_logistic()` |

**Important distinction:** The power simulations use fixed binomial response probabilities to generate trial outcomes and Beta posteriors for analysis. The logistic regression is a separate covariate-adjusted model. Power from the first model is **not automatically** the operating characteristic of the logistic model.

## Installation

Install from the GitHub repository (requires access if it is private):

```r
install.packages("remotes")
remotes::install_github("ASA-BIOP-COMPASS/ST1_2_Bayesian_Design_Ana", upgrade = "never")
library(BayesianDesignAna)
```

Or download the repository ZIP, extract it, and install the folder containing `DESCRIPTION`:

```r
remotes::install_local("~/Downloads/BayesianDesignAna", force = TRUE, upgrade = "never")
library(BayesianDesignAna)
packageVersion("BayesianDesignAna")
```

The **repository name** is `ST1_2_Bayesian_Design_Ana`; the **R package name** used in `library()` is `BayesianDesignAna`. If the repository ZIP extracts to `ST1_2_Bayesian_Design_Ana-main`, use that actual folder path instead. For MAP ESS estimation, also run `install.packages("RBesT")`.

## 1. Set sample sizes and allocation ratios

### `allocation_unit()`, `candidate_totals()`, and `allocate_arms()`

These functions keep the planned sample size in **whole patients** while preserving the **exact** allocation ratio.

```r
ratios <- c(Placebo = 1, Low = 1, Mid = 1, High = 1)

allocation_unit(ratios)                   # 4 patients per allocation block
candidate_totals(ratios, 160, 176)       # 160, 164, 168, 172, 176
allocate_arms(164, ratios)               # 41, 41, 41, 41

# Screen every fifth valid allocation block (5 x 4 = 20 patients):
candidate_totals(ratios, 160, 240, step_blocks = 5)
# 160, 180, 200, 220, 240
```

| Option | Meaning |
| --- | --- |
| `ratios` | Named, positive **integer** allocation weights, with control first; e.g., `c(Placebo=1, Low=2, High=1)`. |
| `min_total`, `max_total` | Inclusive search bounds for the **total randomized** sample size. |
| `step_blocks` | Spacing in allocation blocks, **not individual patients**. Default `1` includes every valid total. |
| `total_n` in `allocate_arms()` | A single total that must be exactly compatible with `ratios`; incompatible totals raise an error. |

For `1:1:1:1`, valid totals are multiples of four (160, 164, 168, ...). For `1:2:2:1`, valid totals are multiples of six (120, 126, 132, ...). The allocation weights are reduced to their simplest integer ratio.

## 2. Estimate power for a fixed design

### `design_power()`

Use this when you **already know the number of patients in each arm** and want the probability of declaring one or more doses successful under assumed response rates.

```r
n <- c(Placebo = 40, Low = 40, High = 40)
rates <- c(Placebo = 0.30, Low = 0.40, High = 0.55)

power <- design_power(
  n = n,
  response_rates = rates,
  borrowed_n = 20,
  historical_rate = 0.30,
  overall = "any",
  multiplicity = "hochberg",
  alpha = 0.025,
  margin = 0,
  n_sim = 2000,
  n_post = 2000,
  seed = 123
)

power$per_arm          # Power and Monte Carlo standard error for each dose
power$overall_power    # Power for the chosen overall success definition
power$any_power        # Probability at least one dose succeeds
power$all_power        # Probability all doses succeed
```

**Key options for power and design functions**

| Option | Meaning | Default / example |
| --- | --- | --- |
| `n` | Named integer patients per arm; control first. | `c(Placebo=40, Low=40, High=40)` |
| `response_rates` | Assumed *true* response probability in each arm; same names and order as `n` or `ratios`. | `c(Placebo=.30, Low=.40, High=.55)` |
| `borrowed_n` | Fixed historical **effective sample size (ESS)** added to the control prior; not randomized patients. | `NULL` (none unless another option is set) |
| `borrowing_prop` | Historical ESS as a fraction of **combined** nominal control information: `ESS / (ESS + concurrent control N)`. | `NULL`; use **either** this or `borrowed_n`, never both. |
| `historical_rate` | Response rate assumed for borrowed historical control information. | `NULL`: uses assumed concurrent control rate in simulations. |
| `overall` | `"any"`: at least one dose succeeds. `"all"`: every dose succeeds. | `"any"` |
| `margin` | Required risk difference: `p_treatment - p_control > margin`. `0` means any positive improvement; `.10` means over 10 percentage points. | `0` |
| `alpha` | One-sided threshold for the **adjusted posterior-tail decision score**. With no adjustment and `alpha=.025`, success requires approximately `Pr(RD > margin | data) >= .975`. | `0.025` |
| `multiplicity` | How posterior-tail scores are adjusted across doses; see below. | `"hochberg"` |
| `n_sim` | Number of simulated trials **per candidate sample size**; more trials reduce Monte Carlo uncertainty. | `500` in `simulate_power()`; `1000` in `design_sample_size()` |
| `n_post` | Posterior Beta draws **per simulated trial**; larger values make the posterior probability approximation more precise, but increase runtime. | `2000` in power simulation |
| `base_control_prior` | Baseline control Beta prior shapes `(alpha, beta)` **before** historical borrowing. | `c(.5,.5)` |
| `treatment_prior` | Beta prior shapes shared by treatment arms. | `c(.5,.5)` |
| `seed` | Reproducible random-number seed. | `123` for design functions |

**Multiplicity choices:** `"none"` (no adjustment), `"bonferroni"` (single-step conservative adjustment), `"holm"` (step-down), `"hochberg"` (step-up), or `"BH"` (Benjamini–Hochberg). These names describe how the package adjusts **Bayesian posterior-tail scores**, not ordinary frequentist p-values. Choosing one **does not by itself guarantee** frequentist family-wise error rate (FWER) or false discovery rate (FDR) control. Evaluate the full decision rule by simulation under null scenarios.

## 3. Find a sample size for target power

### `design_sample_size()`

Use this when you know the **allocation ratio and assumed response rates**, but want to choose the total sample size needed to reach a target power.

```r
ratios <- c(Placebo = 1, Low = 1, Mid = 1, High = 1)
rates <- c(Placebo = .30, Low = .40, Mid = .50, High = .60)
candidates <- candidate_totals(ratios, 120, 240)

ss <- design_sample_size(
  total_n = candidates,
  ratios = ratios,
  response_rates = rates,
  target_power = .80,
  borrowing_prop = .40,
  search = "two_stage",
  initial_step = 5,
  n_sim_initial = 100,
  refine_span = 5,
  n_sim = 1000,
  overall = "any",
  seed = 123,
  plot = TRUE
)

ss$minimum_candidate  # Selected total N and arm allocation, if found
ss$initial_results    # Lower-precision screening
ss$final_results      # Higher-precision results used for selection
ss$results            # Combined results
ss$note               # Any borrowing-related note
plot(ss)              # Replots saved results; no extra simulation
```

**Which search should I choose?**

| Option | Meaning | When to use it |
| --- | --- | --- |
| `search = "grid"` | Runs `n_sim` trials for **every valid** candidate total. | A small candidate range or final confirmation. |
| `search = "two_stage"` | Screens a subset first, then refines around the preliminary target-power crossing. | A broad range where full simulation is costly. |
| `initial_step = 5` | Stage 1 checks every **fifth candidate** in the *full valid candidate grid*. | Faster first-stage screening. |
| `n_sim_initial = 500` | Simulated trials at each Stage 1 candidate. | Small number for screening. |
| `refine_span = 5` | Number of **valid candidate positions** on either side of the preliminary crossing to include in Stage 2 (the code also incorporates an earlier coarse point). | Increase when the power curve is noisy. |
| `n_sim = 10000` | Simulated trials per **Stage 2** candidate (`grid`: per candidate). | Higher precision near crossing. |
| `target_power = .80` | Estimated overall power needed to select a candidate. | Change to your design target. |
| `historical_ess` | Optional historical ESS **reference for a warning only**. It does **not** change borrowing or determine sample size. | Use after an independent MAP ESS analysis. |
| `table_max_cols` | Maximum sample-size columns shown in the plot's allocation/power table. | `12` |
| `plot` | Whether to draw the power curve and allocation table immediately. | `FALSE` |

**Important:** Pass the **full valid grid** to `total_n`. Do not use `step_blocks=5` to remove intermediate candidates before a two-stage search. For a 1:1:1:1 design, a Stage 1 screen might check `160, 180, 200, 220, 240`, but Stage 2 can then check `200, 204, 208, 212, 216, 220`. The selected N is based on **Stage 2 estimated power**, not guaranteed to be the true minimum; confirm nearby candidates with additional simulations.

### `plot_power()` — compare dose-effect scenarios

Each row of `delta_matrix` specifies the **absolute improvement over placebo** for every treatment arm. One row produces one power curve.

```r
deltas <- rbind(
  c(.05, .10, .20),
  c(.10, .20, .30)
)
colnames(deltas) <- c("Low", "Mid", "High")

curves <- plot_power(
  delta_matrix = deltas,
  total_n = candidates,
  ratios = ratios,
  placebo_rate = .30,
  borrowing_prop = .40,
  target_power = .80,
  search = "two_stage",
  n_sim_initial = 100,
  n_sim = 1000,
  table_max_cols = 10,
  plot = TRUE
)
curves$results
```

`placebo_rate` is the assumed placebo response rate; each dose rate equals `placebo_rate + its delta`. `total_n`, borrowing, `overall`, `multiplicity`, priors, simulation settings, and `seed` have the meanings described above. `plot_power()` calls the sample-size search for each scenario.

## 4. Specify historical-control borrowing

These tools serve **different purposes**. Borrowing is never automatic.

### `borrowing_ess()` — convert a borrowing fraction to ESS

```r
borrowing_ess(n_control = 30, borrowing_prop = .40)  # 20
```

Here `n_control=30` is the number of **concurrent randomized controls**. `borrowing_prop=.40` means historical data contribute 40% of combined nominal control information, so `ESS = 30 * .40 / (1 - .40) = 20`. It does **not** mean 40% of 30.

### `borrowing_prior()` — construct a Beta prior

```r
borrowing_prior(historical_rate = .30, ess = 20,
                base_prior = c(.5, .5))
# Beta(6.5, 14.5): .5 + 20*.30, .5 + 20*.70
```

`historical_rate` is the assumed historical response rate; `ess` is borrowed pseudo-patient information; `base_prior` contains the two baseline Beta shape parameters. This **creates prior parameters**; it does not add actual patients to your data.

### `borrowing_limits()` — calculate a planned borrowing cap

```r
borrowing_limits(
  n_control = 30,
  historical_ess = 50,
  max_fraction = .40,
  max_absolute = 25,
  max_ratio = 2
)
# Allowed ESS = min(50, 20, 25, 60) = 20
```

| Option | Meaning |
| --- | --- |
| `historical_ess` | Available historical ESS, such as an RBesT-derived reference. |
| `max_fraction` | Largest historical share of **combined** nominal control information (default `.50`). |
| `max_absolute` | Largest number of borrowed pseudo-patients (default `Inf`, no cap). |
| `max_ratio` | Maximum borrowed ESS divided by concurrent control N (default `Inf`, no cap). |

The result contains `allowed_ess`, `restricted_ess`, and `borrowed_fraction`. This is a **user-specified limit**, not a data-adaptive borrowing model. It is not automatically applied in `design_power()` or `design_sample_size()`; explicitly pass the chosen `borrowed_n`.

### `historical_map_ess()` — estimate MAP prior ESS using RBesT

```r
install.packages(
  "RBesT",
  repos = "https://cloud.r-project.org",
  type = "binary"
)    # Run once if needed
historical <- data.frame(
  study = paste0("Study", 1:6),
  n = c(154, 163, 183, 180, 72, 74),
  r = c(40, 45, 57, 61, 21, 24)
)

map <- historical_map_ess(historical, seed = 123)
map$ess                    # ELIR, moment, and Morita ESS estimates
reference_ess <- min(map$ess)
```

The input columns are `study` (historical study ID), `n` (participants), and `r` (responders). `seed` makes the model fitting reproducible; `...` passes additional arguments to `RBesT::gMAP()`. The function returns the MAP model, its mixture approximation, and ESS estimates using three methods. Using their minimum is an optional conservative **reference choice**, not an automatically validated cap. This function does **not** feed MAP information directly into the design simulator.

## 5. Analyze binary trial outcomes

### `analyze_binary()` — responder counts by arm

Use this after observing **responders and total patients**, without patient-level covariates.

```r
fit_bin <- analyze_binary(
  responses = c(Placebo = 12, Low = 18, High = 25),
  n = c(Placebo = 40, Low = 40, High = 40),
  control_prior = c(.5, .5),
  treatment_prior = c(.5, .5),
  margin = 0, alpha = .025,
  multiplicity = "hochberg",
  n_post = 4000, seed = 123
)
fit_bin
```

`responses` is the number of responders per arm; `n` is the number of observed patients. `control_prior` and `treatment_prior` are Beta shape pairs (not response rates). For explicit borrowing, use `borrowing_prior()` first and supply its `alpha` and `beta` as `control_prior`. Results include posterior mean risk differences (`rd_mean`), 95% **equal-tailed** bounds (`rd_lower`, `rd_upper`), probability of superiority (`prob_superiority`), adjusted posterior-tail scores, and decisions (`success`). `n_post`, `margin`, `alpha`, `multiplicity`, and `seed` have the meanings defined above.


## 6. Analyze patient-level data with covariate adjustment

### `simulate_logistic_data()` — create example patient-level data

```r
n <- c(Placebo=40, Low=40, High=40)
rates <- c(Placebo=.30, Low=.40, High=.55)

dat <- simulate_logistic_data(
  n = n, response_rates = rates,
  strata = list(region = c(North=.6, South=.4)),
  stratum_log_or = list(region = .30),
  seed = 123
)
head(dat)
```

`strata` gives probabilities for each level of a categorical factor; they must sum to 1. `stratum_log_or` gives the log-odds effect for each **nonreference** level (the first listed level is the reference). With more than two levels, supply one effect per nonreference level. The `response_rates` values determine response probabilities at the reference covariate levels, **not necessarily marginal rates across all strata**. Omit `strata` to simulate data without covariates.

### `analyze_logistic()` — fit Bayesian logistic regression

```r
fit <- analyze_logistic(
  data = dat,
  outcome = "y", arm = "arm", control = "Placebo",
  covariates = "region",
  center = "overall",
  control_prior = c(.5, .5),
  interval_method = "both", credible_level = .95,
  sas_compatible_summary = TRUE,
  n_iter = 12000, burnin = 2000,
  thin = 5, n_chains = 2,
  proposal_sd = .15,
  multiplicity = "hochberg", seed = 123
)

fit$summary        # Adjusted risk differences and credible intervals
fit$reference_means
fit$sas_summary    # Optional wide table using SAS-style names
fit$mcmc           # Convergence/ESS summaries
fit$acceptance     # Metropolis proposal acceptance fractions
print(fit)
mcmc_diagnostics(fit$chains)
```

**Logistic regression options**

| Option | Meaning | Default |
| --- | --- | --- |
| `data` | Patient-level data frame with binary outcome (0/1), arm, and optional covariates. | Required |
| `outcome`, `arm`, `control` | Column names for response and assigned treatment, and the control arm's label. | `"y"`, `"arm"`, `"Placebo"` |
| `covariates` | Names of numeric or categorical predictors; factors with *K* levels use *K−1* indicators. | `character()` |
| `center` | `"overall"`: one pooled covariate reference for all arms; `"arm"`: separate arm-specific references; `"none"`: no centering. | `"overall"` |
| `control_prior` | Beta shapes for control response at the reference covariate profile; no historical prior is imposed automatically. | `c(.5,.5)` |
| `intercept_sd`, `covariate_sd` | Standard deviations of Normal(0, SD²) priors for active-arm intercepts and shared covariate coefficients. | `100`, `100` |
| `n_iter`, `burnin`, `thin` | Total iterations per chain, initial iterations discarded, and retention interval. | `12000`, `2000`, `5` |
| `n_chains` | Number of independently initialized MCMC chains. | `2` |
| `proposal_sd` | Random-walk Metropolis proposal step size, not a prior SD. Influences mixing and acceptance. | `.15` |
| `interval_method` | `"equal_tail"`, `"HPD"`, or `"both"`; affects reported interval boundaries only. | `"both"` |
| `credible_level` | Posterior mass included in credible interval. | `.95` |
| `sas_compatible_summary` | Adds `fit$sas_summary` with SAS-style column names and HPD-style bounds; does not make the sampler identical to SAS. | `FALSE` |
| `margin`, `alpha`, `multiplicity`, `seed` | Same threshold, posterior-tail decision rule, adjustment choice, and reproducibility settings described above. | `0`, `.025`, `"hochberg"`, `123` |

**Credible intervals:** Equal-tailed intervals take the lower/upper posterior quantiles (e.g., 2.5% and 97.5% for 95% intervals). The HPD option uses the shortest contiguous interval among retained draws, approximating the HPD interval for a unimodal posterior. With `"both"`, the main `rd_lower`/`rd_upper` columns remain equal-tailed for backward compatibility; separate HPD columns are also returned. The `sas_summary` uses HPD-style limits, and `Pr_Sup` means `Pr(RD > 0 | data)` even if `margin` is nonzero.

**Interpretation:** `center="overall"` evaluates dose and control probabilities at the **same pooled reference covariate profile**. This is a reference-profile risk difference, **not necessarily** a marginal population-average risk difference. `center="arm"` uses different profiles and is not appropriate when the desired estimand is a common-profile contrast. Review MCMC diagnostics and consider sensitivity to priors and starting values.

## 7. Check software and interpret results

For a local source folder:

```r
install.packages("devtools")
devtools::test("~/Downloads/BayesianDesignAna")
devtools::check("~/Downloads/BayesianDesignAna")
```

Use `?design_sample_size`, `?analyze_logistic`, or `help(package="BayesianDesignAna")` for function help and complete signatures.

**Research-use considerations:** This package is development software. Confirm simulation results with larger `n_sim` near the selected N; study global and partial-null operating characteristics of the Bayesian multiplicity decision rule; assess prior–data conflict when borrowing; and validate the custom logistic MCMC sampler before relying on results for confirmatory decisions. A successful R package check does not establish statistical validity.

## References

The historical MAP prior ESS helper builds on **RBesT (R Bayesian Evidence Synthesis Tools)**:

Weber S, Li Y, Seaman JW, Kakizume T, Schmidli H (2021). *Applying Meta-Analytic-Predictive Priors with the R Bayesian Evidence Synthesis Tools.* **Journal of Statistical Software**, 100(19), 1–32. https://doi.org/10.18637/jss.v100.i19

RBesT documentation: https://opensource.nibr.com/RBesT/

License: MIT (see `LICENSE.md`). Report bugs or feature requests at https://github.com/ASA-BIOP-COMPASS/ST1_2_Bayesian_Design_Ana/issues.
