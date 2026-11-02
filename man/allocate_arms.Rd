\name{allocate_arms}
\alias{allocate_arms}
\title{Allocate a total sample size in an exact integer arm ratio}
\usage{allocate_arms(total_n, ratios = c(Control = 1, Dose1 = 1))}
\arguments{
\item{total_n}{Total recruited/randomized subjects; must be a positive integer multiple of the allocation unit.}
\item{ratios}{Named positive integer allocation weights, control first.}
}
\value{Named integer vector, each count exactly satisfying the desired ratio.}
\description{Unlike earlier versions, sample sizes are not rounded to fit incompatible totals. For a 1:1:1:1 ratio, 164 produces 41 subjects per arm and 165 raises an error. Use \code{candidate_totals()} for feasible totals.}
\examples{
allocate_arms(164, c(Placebo=1,Low=1,Mid=1,High=1))
}
