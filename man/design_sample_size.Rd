\name{design_sample_size}
\alias{design_sample_size}
\title{Search exact integer allocations for target Bayesian overall power}
\usage{
design_sample_size(total_n, ratios, response_rates,
  target_power = 0.8, borrowed_n = NULL, borrowing_prop = NULL,
  historical_ess = NULL, plot = FALSE, overall = c("any", "all"),
  multiplicity = "hochberg", seed = 123,
  base_control_prior = c(0.5, 0.5), treatment_prior = c(0.5, 0.5),
  search = c("grid", "two_stage"), n_sim = 1000,
  n_sim_initial = 500, initial_step = 5L, refine_span = NULL,
  table_max_cols = 12L, ...)
}
\arguments{
\item{total_n}{Prespecified randomized totals; totals incompatible with the exact allocation ratio are omitted with a warning. Use \code{candidate_totals()} to generate valid totals.}
\item{ratios}{Named positive integer allocation weights, placebo first.}
\item{response_rates}{Named anticipated response probabilities, same arms and order.}
\item{target_power}{Target overall power, default 0.8.}
\item{borrowed_n}{Fixed historical-control pseudo-subject ESS; mutually exclusive with borrowing_prop.}
\item{borrowing_prop}{Borrowed ESS divided by borrowed ESS plus concurrent placebo count.}
\item{historical_ess}{Independent MAP ESS reference for warnings only; not used in power or sample-size optimization.}
\item{plot}{Plot stored results, without rerunning simulation.}
\item{overall}{Default \code{"any"} means at least one dose succeeds; \code{"all"} means every dose succeeds.}
\item{multiplicity}{Posterior tail-score method: \code{none}, \code{bonferroni}, \code{holm}, \code{hochberg} (default), or \code{BH}. Does not imply frequentist FWER control.}
\item{seed}{Reproducible simulation seed; default 123.}
\item{base_control_prior}{Control Beta prior shapes before borrowing, default c(.5,.5).}
\item{treatment_prior}{Common active-dose Beta prior shapes, default c(.5,.5).}
\item{search}{\code{"grid"} evaluates every feasible candidate with n_sim. \code{"two_stage"} screens a coarse subset and refines around the crossing with n_sim simulations.}
\item{n_sim}{Number of simulated trials per candidate in full search or FINAL stage of two-stage search.}
\item{n_sim_initial}{Number of trials per preliminary coarse-stage candidate.}
\item{initial_step}{Every initial_step-th feasible total is evaluated initially; last total is also evaluated.}
\item{refine_span}{Number of feasible totals on each side of the preliminary crossing/anchor. Defaults to initial_step.}
\item{table_max_cols}{Maximum number of columns in transposed plot table.}
\item{...}{Further arguments passed to design_power(), e.g. n_post, historical_rate, alpha, margin.}
}
\value{An object of class \code{bda_sample_size} with \code{results} (combined estimates, including stage and n_sim), \code{initial_results}, \code{final_results}, \code{minimum_candidate}, allocation information, and optional borrowing reference \code{note}. Final selection uses only high-precision final-stage estimates.}
\description{Power-based integer allocation design with optional resource-saving two-stage Monte Carlo search.}
\details{The two-stage algorithm chooses a preliminary first crossing (or highest-power anchor if none crosses), then evaluates a neighboring window at high precision. It does not guarantee discovery of the mathematical global minimum, especially with Monte Carlo noise or a nonmonotone curve. The historical MAP ESS is not fitted automatically.}
\examples{
r <- c(Placebo=1,Low=1,Mid=1,High=1)
rr <- c(Placebo=.3,Low=.4,Mid=.5,High=.6)
# ss <- design_sample_size(candidate_totals(r,120,200),r,rr,
#   search="two_stage",n_sim_initial=100,n_sim=1000)
}
