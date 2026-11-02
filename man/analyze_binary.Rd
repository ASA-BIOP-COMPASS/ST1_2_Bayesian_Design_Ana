\name{analyze_binary}
\alias{analyze_binary}
\title{Conjugate Bayesian analysis of binary treatment arms}
\usage{
analyze_binary(responses, n, control_prior = c(0.5, 0.5),
  treatment_prior = c(0.5, 0.5), n_post = 4000, margin = 0,
  alpha = 0.025, multiplicity = "hochberg", seed = NULL)
}
\description{Uses independent Beta treatment posteriors and one shared Beta
placebo posterior to estimate risk differences and success probabilities.}
\arguments{
\item{responses}{Named responder counts, placebo first.}
\item{n}{Named patient counts in each arm.}
\item{control_prior}{Two positive placebo Beta shapes; default c(0.5,0.5).}
\item{treatment_prior}{Two positive active-dose Beta shapes; default c(0.5,0.5).}
\item{n_post}{Number of joint posterior draws.}
\item{margin}{Clinical risk-difference threshold.}
\item{alpha}{Posterior-tail threshold.}
\item{multiplicity}{Tail-score adjustment method, default hochberg.}
\item{seed}{Optional random seed.}
}
\value{Data frame with arm-level posterior summaries and decision indicators.}
\details{The default Beta(0.5,0.5) is Jeffreys' reference prior for a
standalone Bernoulli probability. Tail-score adjustment is a Bayesian
operating rule, not automatically a frequentist p-value correction.}
