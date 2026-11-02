\name{simulate_logistic_data}
\alias{simulate_logistic_data}
\title{Simulate patient-level multi-arm binary trial data with stratification factors}
\usage{simulate_logistic_data(n, response_rates, strata = list(), stratum_log_or = list(), seed = 123)}
\description{Generates a patient-level dataset with multiple categorical, possibly multi-level stratification factors and arm-specific baseline response probabilities.}
\value{Data frame containing id, arm, y and stratification factors.}
