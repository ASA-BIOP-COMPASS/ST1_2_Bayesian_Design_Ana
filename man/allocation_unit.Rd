\name{allocation_unit}
\alias{allocation_unit}
\title{Smallest integer allocation block}
\usage{allocation_unit(ratios)}
\arguments{
\item{ratios}{Named positive integer randomization weights, first arm control.}
}
\value{An integer equal to the sum of allocation weights after dividing out their greatest common divisor.}
\description{A trial with an exact allocation ratio must randomize a total that is a positive multiple of this block size. For 1:1:1:1, the block size is 4; for 1:2:2:1, it is 6.}
\examples{
allocation_unit(c(Placebo=1, Low=1, Mid=1, High=1))
}
