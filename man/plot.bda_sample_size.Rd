\name{plot.bda_sample_size}
\alias{plot.bda_sample_size}
\title{Plot stored sample-size search power and arm allocation table}
\usage{\method{plot}{bda_sample_size}(x, ...)}
\arguments{
\item{x}{An object returned by design_sample_size().}
\item{...}{Unused.}
}
\description{Draws the original power estimates and a transposed table of total N, individual arm sizes, borrowed ESS, and overall power. A two-stage search distinguishes preliminary and final estimates and labels the first final total meeting the target power on the x-axis. Rendering does not rerun simulations.}
\value{The input object, invisibly.}
