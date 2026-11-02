\name{candidate_totals}
\alias{candidate_totals}
\title{Generate integer sample-size totals compatible with an exact ratio}
\usage{candidate_totals(ratios, min_total, max_total, step_blocks = 1L)}
\arguments{
\item{ratios}{Named positive integer allocation weights.}
\item{min_total}{Inclusive lower bound for total randomized subjects.}
\item{max_total}{Inclusive upper bound for total randomized subjects.}
\item{step_blocks}{Positive integer. Return every specified number of feasible allocation blocks; default 1.}
}
\value{An integer vector of feasible randomized totals. An empty vector means the range contains no permissible total.}
\description{Only exact-ratio integer allocations are generated.}
\examples{
r <- c(Placebo=1, Low=1, Mid=1, High=1)
candidate_totals(r, 160, 180)
candidate_totals(r, 160, 240, step_blocks=5)
}
