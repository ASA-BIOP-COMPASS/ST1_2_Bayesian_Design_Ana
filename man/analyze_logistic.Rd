\name{analyze_logistic}
\alias{analyze_logistic}
\title{Bayesian multi-arm logistic regression}
\description{Fits Bayesian logistic regression with shared covariate effects and
arm-specific log-odds intercepts using a reference Metropolis sampler.}
\usage{
analyze_logistic(data, outcome = "y", arm = "arm", control = "Placebo",
  covariates = character(), center = "overall",
  control_prior = c(0.5, 0.5), intercept_sd = 100, covariate_sd = 100,
  n_iter = 12000, burnin = 2000, thin = 5, n_chains = 2,
  proposal_sd = 0.15, margin = 0, alpha = 0.025,
  multiplicity = "hochberg", seed = 123,
  interval_method = "both", credible_level = 0.95,
  sas_compatible_summary = FALSE)
}
\arguments{
\item{data}{Patient-level data frame.}
\item{outcome}{Name of the binary outcome column coded as 0 and 1.}
\item{arm}{Name of the treatment arm column.}
\item{control}{Control treatment label.}
\item{covariates}{Zero or more numeric or categorical covariate names.}
\item{center}{Covariate centering method: overall, arm, or none.}
\item{control_prior}{Two positive Beta shapes for the reference placebo probability. Default c(0.5,0.5).}
\item{intercept_sd}{Prior standard deviation for active-dose log-odds intercepts; default 100.}
\item{covariate_sd}{Prior standard deviation for covariate coefficients; default 100.}
\item{n_iter}{Total MCMC iterations per chain, including burn-in.}
\item{burnin}{Number of initial iterations discarded.}
\item{thin}{Retain every thin-th post-burn-in draw.}
\item{n_chains}{Number of independent MCMC chains.}
\item{proposal_sd}{Random-walk proposal standard deviation.}
\item{margin}{Superiority margin on the risk-difference scale.}
\item{alpha}{Posterior tail-score decision threshold.}
\item{multiplicity}{One of none, bonferroni, holm, hochberg, or BH; Hochberg is default.}
\item{seed}{RNG seed.}
\item{interval_method}{Either \code{"both"} (default), \code{"HPD"}, or \code{"equal_tail"}. In \code{"both"} mode the original rd_lower/rd_upper remain equal-tailed for backward compatibility, with extra HPD and equal-tailed columns.}
\item{credible_level}{Credible interval mass strictly between 0 and 1, default 0.95.}
\item{sas_compatible_summary}{Logical; when TRUE, also return a one-row \code{sas_summary} with SAS-style posterior placebo/dose means, risk difference means/HPD bounds, and probabilities of superiority at RD > 0. Dose indices D1, D2, etc. follow the arm order found in the data.}
}
\value{A bbd_logistic list containing treatment effect summaries, posterior draws,
MCMC diagnostics, optional SAS-style one-row summary, the selected interval method and credible level, and prior specifications.}
\details{The default control prior is Beta(0.5,0.5), a reference prior for a
standalone Bernoulli proportion, not a prior learned from historical trials.
Dose and covariate coefficients receive diffuse proper Normal(0,100^2)
priors by default. An informative placebo prior may be provided through
\code{control_prior} when justified by suitable evidence. The very diffuse Normal
priors are not the joint objective prior for this full logistic model.
HPD intervals approximate shortest contiguous intervals of the retained draws and may differ from SAS for finite MCMC sampling or multimodal distributions. Convergence and frequentist operating characteristics require validation.}
