\name{design_power}
\alias{design_power}
\title{Calculate Bayesian power for a multi-arm binary trial}
\usage{
design_power(n, response_rates, borrowed_n = NULL, borrowing_prop = NULL,
  historical_rate = NULL, overall = c("any", "all"),
  multiplicity = "hochberg", seed = 123,
  base_control_prior = c(0.5, 0.5), treatment_prior = c(0.5, 0.5), ...)
}
\description{Computes power from one concurrent control and flexible active-dose arms.}
\arguments{
\item{n}{Named recruited arm sample sizes, placebo first.}
\item{response_rates}{Named anticipated response probabilities.}
\item{borrowed_n}{Optional fixed historical borrowing ESS; mutually exclusive with borrowing_prop.}
\item{borrowing_prop}{Optional historical fraction of combined placebo information.}
\item{historical_rate}{Historical placebo response probability when borrowing.}
\item{overall}{any (default) or all active doses must succeed.}
\item{multiplicity}{Posterior-tail multiplicity method; default hochberg.}
\item{seed}{Simulation seed.}
\item{base_control_prior}{Baseline Beta control shapes; default c(0.5,0.5).}
\item{treatment_prior}{Beta shapes for each dose; default c(0.5,0.5).}
\item{...}{Simulation controls passed to simulate_power(), such as n_sim and n_post.}
}
\value{List with per-arm power, overall_power, borrowed_ess, and scenario.}
\details{Borrowing is off by default. If requested, borrowed placebo
information is added to the baseline control Beta shapes. Hochberg adjusts
Bayesian posterior-tail decision scores, not conventional frequentist p-values.}
