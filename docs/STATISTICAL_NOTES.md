# Statistical notes

Design: binomial sampling with independent Beta posterior treatment probabilities and one shared Beta control posterior, shaped by a fixed borrowing proportion. This is distinct from the historical-study MAP model and from covariate-adjusted Bayesian logistic regression.

Overall any-dose power is the probability at least one dose is declared successful; all-dose power requires every dose declared successful. Dose-specific power is marginal across simulation replicates.

Posterior tail scores after adjustment are decision conventions; their use in p.adjust does not guarantee frequentist FWER control. Simulate null and alternative configurations to calibrate thresholds. Borrowing proportion is nominal borrowed ESS / (concurrent control n + nominal borrowed ESS), with a base Jeffreys Beta(0.5,0.5) prior added. Historical MAP ESS can be calculated separately but does not directly validate a borrowing proportion.

`plot_power` accepts a delta matrix: each row varies the entire vector of active-dose effects and produces one curve. The design engine uses fixed binomial response rates across trials; a Beta-binomial data generator would specify a different operating-characteristic model. The logistic sampler is a custom Metropolis implementation and needs convergence/benchmark validation before scientific reliance.

## Prior policy

The default unborrowed control and active-dose priors in binary design are
Beta(0.5, 0.5). The logistic analysis uses Beta(0.5,0.5) for the reference
placebo response probability and independent, very diffuse proper
Normal(0, 100^2) priors for dose log-odds and shared covariate effects.
Prior sensitivity is essential: Jeffreys' Beta prior is invariant for a
standalone Bernoulli proportion but does not become the objective Jeffreys
prior for the entire adjusted logistic model. A diffuse Normal can be weak
informationally yet cause mixing issues or instabilities with separation.
Borrowing is off by default and must be explicitly requested. Historical
MAP model hyperpriors are separate and have not been changed.

## Posterior credible intervals

The `analyze_logistic()` function supports equal-tailed quantile intervals and
shortest contiguous posterior-draw intervals (HPD approximations). When
`interval_method="both"`, `rd_lower` and `rd_upper` are equal-tailed for backward
compatibility, and both named pairs are returned. HPD and equal-tail differ for
skewed posteriors; neither choice changes posterior means or superiority
probabilities. The optional `sas_summary` uses HPD-style interval limits,
and is intended as a comparison aid rather than a claim of numerical equivalence
with SAS PROC MCMC.
