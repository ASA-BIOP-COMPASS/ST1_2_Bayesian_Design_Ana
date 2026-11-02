\name{find_sample_size}
\alias{find_sample_size}
\title{Legacy sample-size simulation search over exact-ratio candidates}
\usage{find_sample_size(total_n, ratios, response_rates, target_power = 0.8, target = "any", ...)}
\arguments{
\item{total_n}{Integer total randomized candidates; infeasible totals are skipped, since allocations must exactly preserve the ratio.}
\item{ratios}{Named positive integer allocation weights with control first.}
\item{response_rates}{Named treatment and control response probabilities aligned with ratios.}
\item{target_power}{Minimum desired power, default 0.8.}
\item{target}{'any' for at least one success, 'all' for all active arms, or one active arm name.}
\item{...}{Additional arguments forwarded to simulate_power(), including borrow_ess, historical_rate, multiplicity, alpha, n_sim, n_post, and seed.}
}
\value{List containing results (candidate totals, power estimates, per-arm allocations, success flags) and minimum_candidate.}
\description{Legacy helper for a simple exact-ratio candidate search. For borrowing proportion, historical ESS reference warnings, transposed plots or two-stage searching, prefer design_sample_size().}
\details{Historical MAP ESS is not fitted automatically. Bayesian posterior-tail multiplicity adjustments have not been shown to control nominal frequentist error.}
