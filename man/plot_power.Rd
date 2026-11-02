\name{plot_power}
\alias{plot_power}
\title{Compare Bayesian overall power across dose-effect scenarios}
\usage{
plot_power(delta_matrix, total_n, ratios, placebo_rate,
  borrowed_n = NULL, borrowing_prop = NULL, historical_ess = NULL,
  plot = TRUE, overall = c("any", "all"), multiplicity = "hochberg",
  seed = 123, base_control_prior = c(0.5, 0.5),
  treatment_prior = c(0.5, 0.5), table_max_cols = 12L, ...)
}
\arguments{
\item{delta_matrix}{Numeric matrix: one row per scenario and one column per active treatment arm. Values are risk differences versus placebo.}
\item{total_n}{Prespecified candidate randomized totals. Only feasible exact-ratio integer allocations are retained.}
\item{ratios}{Named integer allocation weights, placebo first.}
\item{placebo_rate}{Assumed placebo response probability.}
\item{borrowed_n}{Fixed borrowed placebo ESS; alternative to borrowing_prop.}
\item{borrowing_prop}{Borrowed fraction of total control information.}
\item{historical_ess}{Historical MAP ESS reference for warning only.}
\item{plot}{Whether to draw a curve for each row of the delta matrix.}
\item{overall}{Default 'any' or 'all' dose success rule.}
\item{multiplicity}{Posterior tail-score method, default 'hochberg'.}
\item{seed}{Simulation random seed.}
\item{base_control_prior}{Baseline placebo Beta shapes before borrowing.}
\item{treatment_prior}{Treatment Beta shapes.}
\item{table_max_cols}{Maximum transposed table width, default 12 scenario columns.}
\item{...}{Passed to design_sample_size(), including search, n_sim, n_sim_initial, initial_step, refine_span, target_power, n_post, margin, alpha and historical_rate.}
}
\description{Each scenario calls design_sample_size(), then plots retained simulation results without rerunning work. Scenario-specific two-stage searches may produce different evaluated grids; missing table values are explicitly labeled instead of interpolated.}
\value{List with results (including scenario labels), scenarios (per-scenario search objects), and legend labels, returned invisibly when plot=TRUE.}
