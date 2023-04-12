#' @title Zero-Inflated Negative Binomial Regression
#' @description Conduct zero-inflated negative binomial regression for data. Estimate paramters and their covariance using EM algorithm.
#' @param X A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates. Number of observations \eqn{n} must be greater than number of covariates \eqn{p}.
#' @param y An \eqn{n}-dimensional response vector.
#' @param tol (\strong{optional}) The tolerance level in the iterative estimation procedure. The defalut value is 1e-4.
#' @param maxIter (\strong{optional}) Maximum number of iterations. Default is 1000.
#' @param initial (\strong{optional}) A list specifies the initial point of the algorithm \code{list(beta=..., gamma=..., theta=...)}. Default is NULL.
#' @return A list including the following terms will be returned:
#' \describe{
#' \item{\code{beta}}{The coefficients for the Negative Binomial Model.}
#' \item{\code{gamma}}{The coefficients for the Zero Inflated Model.}
#' \item{\code{theta}}{The reciprocal of the size parameter for the Negative Binomial Distribution.}
#' \item{\code{vcov}}{The covariance matrix for all estimated paramters.}
#' \item{\code{loglikelihood}}{The log-likelihood of the estimated model.}
#' }
#' @examples
#' n <- 100
#' X <- matrix(c(rep(1, 50), rbinom(50, 1, 0.5)), 50, 2)
#' y <- sample(c(rnbinom(20, 10, 0.5), rep(0, 30)))
#' fit.mod <- fit.zinb(X, y)
#' fit.mod$beta
#' @export
#' @useDynLib ZINB
fit.zinb <- function(X, y, tol = 1e-4, maxIter = 1000, initial = NULL) {
  return(par_EM(X, y, tol, maxIter, initial))
}

#' @title Likelihood Ratio Test for Zero-Inflated Negative Binomial Regression (Single Response)
#' @description Testing the intercept model against the model fitted by the design matrix. It is recommended to include the intercept in X.
#' @param X A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates. Number of observations \eqn{n} must be greater than number of covariates \eqn{p}.
#' @param y An \eqn{n}-dimensional response vector.
#' @param ... (\strong{optional}) Other optional parameters for \code{fit.zinb}.
#' @return A list including the following terms will be returned:
#' \describe{
#' \item{\code{par}}{Parameters of the fitted model.}
#' \item{\code{Chisq}}{Calculated Chisq statistics.}
#' \item{\code{df}}{The degree of freedom of this test.}
#' \item{\code{p.value}}{Computed p-values.}
#' }
#' @examples
#' n <- 100
#' X <- matrix(c(rep(1, 50), rbinom(50, 1, 0.5)), 50, 2)
#' y <- sample(c(rnbinom(20, 10, 0.5), rep(0, 30)))
#' result <- LRT1D(X, y)
#' result$p.value
#' @export
LRT1D <- function(X, y, ...){
  return(LRT1G(X, y, ...))
}

#' @title Likelihood Ratio Test for Zero-Inflated Negative Binomial Regression (Multiple Responses)
#' @description Testing the intercept model against the model fitted by the design matrix. It is recommended to include the intercept in X.
#' @param X A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates. Number of observations \eqn{n} must be greater than number of covariates \eqn{p}.
#' @param Y An \eqn{n} by \eqn{p} response matrix.
#' @param ... (\strong{optional}) Other optional parameters for \code{fit.zinb}.
#' @return A list including the following terms will be returned:
#' \describe{
#' \item{\code{Chisq}}{A vector of calculated Chisq statistics.}
#' \item{\code{df}}{The degree of freedom of this test.}
#' \item{\code{p.value}}{A vector of computed p-values.}
#' }
#' @examples
#' n <- 100
#' X <- matrix(c(rep(1, 50), rbinom(50, 1, 0.5)), 50, 2)
#' y <- sample(c(rnbinom(20, 10, 0.5), rep(0, 30)))
#' Y <- matrix(c(y, sample(y)), 50, 2)
#' result <- LRTnD(X, Y)
#' result$p.value
#' @export
LRTnD <- function(X, y, ...){
  return(LRTnG(X, y, ...))
}