% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/foldedSIR-old.R
\name{foldedSIR.old}
\alias{foldedSIR.old}
\title{Old Implementation of Folded Sliced Inverse Regression}
\usage{
foldedSIR.old(y, x, p, nslice, d, maxiter = 200, tol = 1e-08)
}
\arguments{
\item{y}{Numeric vector of response variables.}

\item{x}{Numeric matrix of predictors.}

\item{p}{Numeric vector specifying the tensor dimensions.}

\item{nslice}{Number of slices for SIR.}

\item{d}{Dimensionality reduction targets.}

\item{maxiter}{Maximum iterations for convergence.}

\item{tol}{Tolerance for convergence.}
}
\value{
A list containing model fit information, including parameters and convergence details.
}
\description{
This function is an older version of the \code{foldedSIR} method for dimensionality reduction,
supporting only specific dimensions (\code{d = c(1, 1)}) for bidimensional tensorial data.
}
\examples{
# Example usage would be similar to `foldedSIR`, but with `d` restricted to `c(1, 1)`.
}
