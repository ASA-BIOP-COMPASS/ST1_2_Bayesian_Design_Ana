\name{borrowing_prior}
\alias{borrowing_prior}
\title{Create a Beta prior with a fixed amount of historical borrowing}
\description{Create a Beta prior with a fixed amount of historical borrowing. See the README and docs/STATISTICAL_NOTES.md in the source repository for methodological assumptions and limitations.}
\usage{
borrowing_prior(historical_rate, ess = 0, base_prior = c(0.5, 0.5))
}
\arguments{
\item{...}{See function signature and README for details.}
}
\value{A list with Beta shape parameters and incremental borrowing ESS.}
\details{This software is a research preview. Bayesian posterior-tail adjustment is not guaranteed frequentist multiplicity control. Calibrate and validate before use in clinical studies.}
