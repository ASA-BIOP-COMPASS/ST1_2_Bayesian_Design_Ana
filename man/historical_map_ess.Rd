\name{historical_map_ess}
\alias{historical_map_ess}
\title{Estimate MAP prior effective sample size using RBesT}
\description{Estimate MAP prior effective sample size using RBesT. See the README and docs/STATISTICAL_NOTES.md in the source repository for methodological assumptions and limitations.}
\usage{
historical_map_ess(data, seed = 123, ...)
}
\arguments{
\item{...}{See function signature and README for details.}
}
\value{List with MAP fit, mixture approximation and elir/moment/morita ESS.}
\details{This software is a research preview. Bayesian posterior-tail adjustment is not guaranteed frequentist multiplicity control. Calibrate and validate before use in clinical studies.}
