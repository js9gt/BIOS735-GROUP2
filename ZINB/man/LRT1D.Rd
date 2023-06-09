% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zinb.R
\name{LRT1D}
\alias{LRT1D}
\title{Likelihood Ratio Test for Zero-Inflated Negative Binomial Regression (Single Response)}
\usage{
LRT1D(X, y, add_intercept = TRUE, ...)
}
\arguments{
\item{X}{A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates. Number of observations \eqn{n} must be greater than number of covariates \eqn{p}.}

\item{y}{An \eqn{n}-dimensional response vector.}

\item{add_intercept}{(\strong{optional}) Whether an intercept term should be added to the design matrix X. If so, the first element of \code{beta} and \code{gamma} would be the intercept. Default is TRUE.}

\item{...}{(\strong{optional}) Other optional parameters for \code{fit.zinb}.}
}
\value{
A list including the following terms will be returned:
\describe{
\item{\code{par}}{Parameters of the fitted model.}
\item{\code{Chisq}}{Calculated Chisq statistics.}
\item{\code{df}}{The degree of freedom of this test.}
\item{\code{p.value}}{Computed p-values.}
}
}
\description{
Testing the intercept model against the model fitted by the design matrix. It is recommended to include the intercept in X or set \code{add_intercept} as TRUE.
}
\examples{
n <- 100
X <- matrix(c(rep(1, 50), rbinom(50, 1, 0.5)), 50, 2)
y <- sample(c(rnbinom(20, 10, 0.5), rep(0, 30)))
result <- LRT1D(X, y, add_intercept=FALSE)
result$p.value
}
