#' @importFrom Rcpp evalCpp
#' @import RcppArmadillo
#' @title Zero-Inflated Negative Binomial Regression
#' @description Conduct zero-inflated negative binomial regression for data. Estimate paramters and their covariance using EM algorithm.
#' @param X A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates. Number of observations \eqn{n} must be greater than number of covariates \eqn{p}.
#' @param y An \eqn{n}-dimensional response vector.
#' @param add_intercept (\strong{optional}) Whether an intercept term should be added to the design matrix X. If so, the first element of \code{beta} and \code{gamma} would be the intercept. Default is TRUE.
#' @param tol (\strong{optional}) The tolerance level in the iterative estimation procedure. The defalut value is 1e-3.
#' @param inflation (\strong{optional}) Whether the inflation model is used. Default is TRUE.
#' @param maxIter (\strong{optional}) Maximum number of iterations. Default is 5000.
#' @param initial (\strong{optional}) A list specifies the initial point of the algorithm \code{list(beta=..., gamma=..., theta=...)}. Default is NULL.
#' @return A list including the following terms will be returned:
#' \describe{
#' \item{\code{infaltion}}{Logical variable to indicate whether the model is inflated. When there is no zero count in y, this should be 0.}
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
#' fit.mod <- fit.zinb(X, y, add_intercept = FALSE)
#' fit.mod$beta
#' @export
#' @useDynLib ZINB
fit.zinb <- function(X, y, add_intercept = TRUE, tol = 1e-3, maxIter = 5000, inflation = TRUE, initial = NULL) {
  # check the data structure, if possible make it compliable to par_EM
  X <- as.matrix(X)
  y <- as.vector(y)
  if(add_intercept){
    X <- cbind(1, X)
  }

  result <- par_EM(X, y, tol, maxIter, inflation, initial)

  # make the armadillo more consistent with R standard data structures
  result$beta <- as.vector(result$beta)
  result$gamma <- as.vector(result$gamma)
  return(result)
}

#' @title Likelihood Ratio Test for Zero-Inflated Negative Binomial Regression (Single Response)
#' @description Testing the intercept model against the model fitted by the design matrix. It is recommended to include the intercept in X or set \code{add_intercept} as TRUE.
#' @param X A \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates. Number of observations \eqn{n} must be greater than number of covariates \eqn{p}.
#' @param y An \eqn{n}-dimensional response vector.
#' @param add_intercept (\strong{optional}) Whether an intercept term should be added to the design matrix X. If so, the first element of \code{beta} and \code{gamma} would be the intercept. Default is TRUE.
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
#' result <- LRT1D(X, y, add_intercept=FALSE)
#' result$p.value
#' @export
LRT1D <- function(X, y, add_intercept=TRUE, ...){
  # check the data structure, if possible make it compliable to par_EM
  X <- as.matrix(X)
  y <- as.vector(y)
  if(add_intercept){
    X <- cbind(1, X)
  }

  result <- LRT1G(X, y, ...)

  # correct the data structure of some outputs
  result$par$beta <- as.vector(result$par$beta)
  result$par$gamma <- as.vector(result$par$gamma)
  return(result)
}

#' @importFrom stats p.adjust
#' @title Likelihood Ratio Test for Zero-Inflated Negative Binomial Regression (Multiple Responses)
#' @description Testing the intercept model against the model fitted by the design matrix. It is recommended to include the intercept in X or set \code{add_intercept} as TRUE.
#' @param X An \eqn{n} by \eqn{p} design matrix. Each row is a vector of observation with \eqn{p} covariates. Number of observations \eqn{n} must be greater than number of covariates \eqn{p}.
#' @param Y An \eqn{n} by \eqn{q} response matrix.
#' @param inflation (\strong{optional}) A \eqn{q}-dimensional logical vector indicating whether the inflated model should be used for each test. If set to NULL, use the inflated model for all tests. Default is NULL.
#' @param add_intercept (\strong{optional}) Whether an intercept term should be added to the design matrix X. If so, the first element of \code{beta} and \code{gamma} would be the intercept. Default is TRUE.
#' @param ... (\strong{optional}) Other optional parameters for \code{fit.zinb}.
#' @return A list including the following terms will be returned:
#' \describe{
#' \item{\code{Chisq}}{A vector of calculated Chisq statistics.}
#' \item{\code{df}}{The degree of freedom of the tests.}
#' \item{\code{p.value}}{A vector of computed p-values.}
#' \item{\code{q.value}}{A vector of adjusted p-values using the Benjamini-Hochberg (BH) correction.}
#' }
#' @examples
#' n <- 100
#' X <- matrix(c(rep(1, 50), rbinom(50, 1, 0.5)), 50, 2)
#' y <- sample(c(rnbinom(20, 10, 0.5), rep(0, 30)))
#' Y <- matrix(c(y, sample(y)), 50, 2)
#' result <- LRTnD(X, Y, add_intercept=FALSE)
#' result$p.value
#' @export
LRTnD <- function(X, Y, inflation=NULL, add_intercept=TRUE, ...){
  # check the data structure, if possible make it compliable
  X <- as.matrix(X)
  Y <- as.matrix(Y)

  if(is.null(inflation)){
    inflation <- rep(1, dim(Y)[2])
  }
  else{
    inflation <- as.vector(inflation)
  }
  if(add_intercept){
    X <- cbind(1, X)
  }

  result <- LRTnG(X, Y, inflation, ...)

  # correct the data structure
  result$Chisq <- as.vector(result$Chisq)
  result$df <- as.vector(result$df)
  result$p.value <- as.vector(result$p.value)
  result$q.value <- p.adjust(result$p.value, method="BH")
  return(result)
}
