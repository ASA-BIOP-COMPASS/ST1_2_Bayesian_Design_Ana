\name{mcmc_diagnostics}
\alias{mcmc_diagnostics}
\title{Calculate basic split R-hat and MCMC effective sample size}
\description{Calculate basic split R-hat and MCMC effective sample size. See the README and docs/STATISTICAL_NOTES.md in the source repository for methodological assumptions and limitations.}
\usage{
mcmc_diagnostics(chains)
}
\arguments{
\item{...}{See function signature and README for details.}
}
\value{Data frame with approximate split R-hat and MCMC ESS by parameter.}
\details{This software is a research preview. Bayesian posterior-tail adjustment is not guaranteed frequentist multiplicity control. Calibrate and validate before use in clinical studies.}
