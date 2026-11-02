\name{borrowing_limits}
\alias{borrowing_limits}
\title{Calculate capped borrowed effective sample size}
\description{Calculate capped borrowed effective sample size. See the README and docs/STATISTICAL_NOTES.md in the source repository for methodological assumptions and limitations.}
\usage{
borrowing_limits(n_control, historical_ess, max_fraction = 0.5, max_absolute = Inf, max_ratio = Inf)
}
\arguments{
\item{...}{See function signature and README for details.}
}
\value{One-row data frame describing available and allowable borrowing ESS.}
\details{This software is a research preview. Bayesian posterior-tail adjustment is not guaranteed frequentist multiplicity control. Calibrate and validate before use in clinical studies.}
