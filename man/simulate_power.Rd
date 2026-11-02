\name{simulate_power}
\alias{simulate_power}
\title{Simulate multi-arm Bayesian binary trial operating characteristics}
\description{Simulates trials with a shared placebo posterior and independent
active-dose Beta priors.}
\usage{
simulate_power(n, response_rates, borrow_ess = 0, historical_rate = NULL,
  n_sim = 500, n_post = 2000, multiplicity = "hochberg",
  alpha = 0.025, margin = 0, seed = 123,
  base_control_prior = c(0.5, 0.5), treatment_prior = c(0.5, 0.5))
}
\arguments{
\item{n}{Named sample sizes, with placebo first.}
\item{response_rates}{Named true arm response probabilities.}
\item{borrow_ess}{Historical placebo borrowing pseudo-sample size; default zero.}
\item{historical_rate}{Prior mean for historical placebo borrowing; defaults to concurrent placebo rate when borrowing.}
\item{n_sim}{Number of trial simulation replicates.}
\item{n_post}{Posterior Monte Carlo draws per trial.}
\item{multiplicity}{Decision adjustment method; default Hochberg.}
\item{alpha}{One-sided posterior-tail decision threshold.}
\item{margin}{Superiority margin on risk-difference scale.}
\item{seed}{RNG seed.}
\item{base_control_prior}{Baseline placebo Beta shapes before borrowing; default c(0.5,0.5).}
\item{treatment_prior}{Beta shapes for each treatment arm; default c(0.5,0.5).}
}
\value{List with per-arm, any-dose, and all-dose power and scenario definitions.}
\details{No borrowing is automatic. If borrow_ess is positive, historical
pseudo-information is added to the selected base_control_prior. Adjusting
Bayesian tail scores does not ensure conventional familywise error control.}
