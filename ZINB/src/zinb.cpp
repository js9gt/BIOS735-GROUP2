#define ARMA_WARN_LEVEL 1
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
arma::vec c_pi(const arma::mat& X, const arma::vec& gammak) {
  return exp(X * gammak) / (1 + exp(X * gammak));
}

// [[Rcpp::export]]
arma::vec c_mu(const arma::mat& X, const arma::vec& betak) {
  return exp(X * betak);
}

// [[Rcpp::export]]
arma::vec Qik(const arma::vec& pi_k, const arma::vec& mu_k, const arma::vec& y_k, double thetak) {
  arma::vec indicator = arma::conv_to<arma::vec>::from(y_k == 0);
  return pi_k % indicator / (pi_k + (1 - pi_k) / pow(1 + thetak * mu_k, 1 / thetak));
}

// Score function: Compute the score of the function
// [[Rcpp::export]]
arma::vec score(const arma::mat& X, const arma::vec& y_k, const arma::vec& z_k, const arma::vec& mu_k, const arma::vec& pi_k, double thetak) {
    int n = X.n_rows;
    int p = X.n_cols;

    arma::mat V = arma::zeros<arma::mat>(n, n);
    V.diag() = 1 / (1 + thetak * mu_k);
    arma::vec ey = (1 - z_k) % (y_k - mu_k);
    arma::vec ez = z_k - pi_k;

    arma::vec result(2 * p);
    result.subvec(0, p - 1) = X.t() * V * ey;
    result.subvec(p, 2 * p - 1) = X.t() * ez;

    return result;
}

// Compute the Fisher information
// [[Rcpp::export]]
arma::mat In(const arma::mat& X, const arma::vec& y_k, const arma::vec& z_k, const arma::vec& mu_k, const arma::vec& pi_k, const arma::vec& gammak, double thetak) {
  int n = X.n_rows;
  int p = X.n_cols;

  arma::mat Py = arma::zeros<arma::mat>(n, n);
  arma::mat Pz = arma::zeros<arma::mat>(n, n);
  Py.diag() = (1 - z_k) % mu_k / (1 + thetak * mu_k);
  Pz.diag() = pi_k / (1 + exp(X * gammak));

  arma::mat I_mat = arma::zeros<arma::mat>(2 * p, 2 * p);
  I_mat.submat(0, 0, p - 1, p - 1) = X.t() * Py * X;
  I_mat.submat(p, p, 2 * p - 1, 2 * p - 1) = X.t() * Pz * X;

  return I_mat;
}

// Moment equation of theta
// [[Rcpp::export]]
double M_estimator(double thetak, const arma::vec& yik, const arma::vec& muik, const arma::vec& zik, int p) {
  return arma::sum((1 - zik) % arma::pow(yik - muik, 2) / muik / (1 + thetak * muik) - (1 - zik)) + p;
}

// Derivatives of the moment derivation
// [[Rcpp::export]]
double M_deriv(double thetak, const arma::vec& yik, const arma::vec& muik, const arma::vec& zik) {
  return - arma::sum((1 - zik) % arma::pow(yik - muik, 2) / arma::pow(1 + thetak * muik, 2));
}

// update theta
// [[Rcpp::export]]
double theta_update(double thetak, int numIter, const arma::vec& yik, const arma::vec& muik, const arma::vec& zik, int p) {
  // if the newton method returns something smaller than 0, this indicates that the model is not overdispersed
  // hence we use very small dispersion parameter, the result is close to Poisson regression
  if (M_estimator(0, yik, muik, zik, p) <= 0) {
    thetak = 1e-10;
  } else {
    for (int i = 0; i < numIter; ++i) {
      thetak -= M_estimator(thetak, yik, muik, zik, p) / M_deriv(thetak, yik, muik, zik);
    }
  }
  if (thetak < 0) {
    thetak = 1e-10;
  }
  return thetak;
}

// The density of zero inflated negative binomial distribution
// [[Rcpp::export]]
arma::vec dznbinom(const arma::vec& x, double size, const arma::vec& prob, arma::vec& infl) {
  int n = x.n_elem;
  arma::vec f(n);
  for(int i = 0; i < n; i++) {
    f(i) = R::dnbinom(x(i), size, prob(i), 0);
  }
  arma::vec result = pow((infl + (1 - infl) % pow(prob, size)), arma::conv_to<arma::vec>::from(x == 0)) %
    pow((1 - infl) % f, arma::conv_to<arma::vec>::from(x > 0));
  return result;
}

// EM algorithm
// [[Rcpp::export]]
List par_EM(const arma::mat& X, const arma::vec& y_k, double tol = 1e-3, int maxIter = 5000, int inflation = true, Nullable<List> initial = R_NilValue) {
  if (X.n_rows != y_k.n_elem) {
    stop("The dimension of X and y does not match!");
  }

  if (min(y_k) < 0){
    stop("Invalid input of y. Must be positive or zero.");
  }
  int n = X.n_rows;
  int p = X.n_cols;

  arma::vec betak, gammak;
  double thetak, loglik_new, loglik_old;
  loglik_old = INFINITY;

  // model initialization
  if (initial.isNull()) {
    betak = gammak = arma::zeros<arma::vec>(p);
    thetak = 1;
  } else {
    List initial_list(initial);
    betak = as<arma::vec>(initial_list["beta"]);
    gammak = as<arma::vec>(initial_list["gamma"]);
    thetak = as<double>(initial_list["theta"]);
  }

  // determine whether inflated model should be fitted
  if(min(y_k) != 0) {
    inflation = false;
  }

  int iter = 0;
  double eps = INFINITY;
  arma::mat I_mat, identity;

  arma::vec pi_k(n), mu_k(n), z_k(n);

  // zero inflated model
  if(inflation){
    identity = arma::eye(2 * p, 2 * p);
    I_mat = arma::eye(2 * p, 2 * p);
    arma::vec s(2 * p), par_old(2 * p), par_new(2 * p);
    while (iter < maxIter && eps > tol) {
      // ---- E step ----
      pi_k = c_pi(X, gammak);
      mu_k = c_mu(X, betak);
      z_k = Qik(pi_k, mu_k, y_k, thetak);
      // theta is updated here
      thetak = theta_update(thetak, 15, y_k, mu_k, z_k, p);

      // ---- Evaluating the Model ----
      arma::vec p_0 = 1 / (1 + thetak * mu_k);
      arma::vec logprob = log(dznbinom(y_k, 1 / thetak, p_0, pi_k));
      // remove NA and Infs
      loglik_new = accu(logprob.elem(find_finite(logprob)));
      eps = std::abs(loglik_new - loglik_old);
      if (std::isnan(eps)) {
        eps = INFINITY;
      }
      loglik_old = loglik_new;      

      // ---- M step ----
      s = score(X, y_k, z_k, mu_k, pi_k, thetak);
      I_mat = In(X, y_k, z_k, mu_k, pi_k, gammak, thetak);

      par_old = arma::join_vert(betak, gammak);
      par_new = par_old + arma::solve(I_mat, s, arma::solve_opts::fast + arma::solve_opts::allow_ugly);

      betak = par_new.subvec(0, p - 1);
      gammak = par_new.subvec(p, 2 * p - 1);

      iter = iter + 1;
    }
    // compute the final loglikelihood
    pi_k = c_pi(X, gammak);
    mu_k = c_mu(X, betak);
    arma::vec p_0 = 1 / (1 + thetak * mu_k);
    arma::vec logprob = log(dznbinom(y_k, 1 / thetak, p_0, pi_k));
    // remove NA and Infs
    loglik_new = accu(logprob.elem(find_finite(logprob)));
  }
  else{
    // this is the given information for the non-inflated model
    identity = arma::eye(p, p);
    I_mat = arma::eye(p, p);
    pi_k = arma::zeros(n);
    z_k = arma::zeros(n);
    arma::vec s(p);
    arma::uvec index = arma::regspace<arma::uvec>(0, p - 1);

    // fit a negative binomial regression
    while (iter < maxIter && eps > tol) {
      // compute theta and model evaluating
      mu_k = c_mu(X, betak);
      thetak = theta_update(thetak, 15, y_k, mu_k, z_k, p);
      arma::vec p_0 = 1 / (1 + thetak * mu_k);
      arma::vec logprob = log(dznbinom(y_k, 1 / thetak, p_0, pi_k));
      loglik_new = accu(logprob.elem(find_finite(logprob)));
      eps = std::abs(loglik_new - loglik_old);
      if (std::isnan(eps)) {
        eps = INFINITY;
      }
      loglik_old = loglik_new;  

      // one step Newton-Raphson method
      s = score(X, y_k, z_k, mu_k, pi_k, thetak)(index);
      I_mat = In(X, y_k, z_k, mu_k, pi_k, gammak, thetak)(index, index);
      betak = betak + arma::solve(I_mat, s, arma::solve_opts::fast + arma::solve_opts::allow_ugly);
      iter = iter + 1;
    }
    // compute the final loglikelihood
    mu_k = c_mu(X, betak);
    arma::vec p_0 = 1 / (1 + thetak * mu_k);
    arma::vec logprob = log(dznbinom(y_k, 1 / thetak, p_0, pi_k));
    // remove NA and Infs
    loglik_new = accu(logprob.elem(find_finite(logprob)));
  }

  if (eps > tol) {
    Rcpp::warning("The algorithm completes without convergence.");
  }
  return List::create(
    _["inflation"] = inflation,
    _["beta"] = betak,
    _["gamma"] = gammak,
    _["theta"] = thetak,
    _["vcov"] = arma::solve(I_mat, identity, arma::solve_opts::allow_ugly),
    _["loglikelihood"] = loglik_new,
    _["iterations"] = iter,
    _["error"] = eps
  );
}

// Test Significance for single species
// [[Rcpp::export]]
List LRT1G(const arma::mat& X, const arma::vec& y_k, double tol = 1e-3, int maxIter = 5000, int inflation = true, Nullable<List> initial = R_NilValue) {
  if (X.n_rows != y_k.n_elem) {
    stop("The dimension of X and y does not match!");
  }
  int n = X.n_rows;
  int p = X.n_cols;

  int df;

  // alternative hypothesis
  List par_alt = par_EM(X, y_k, tol, maxIter, inflation, initial);
  double l_alt = as<double>(par_alt["loglikelihood"]);

  // null hypothesis
  arma::mat X_0 = arma::ones<arma::mat>(n, 1);
  List par_0 = par_EM(X_0, y_k, tol, maxIter, inflation, initial);
  double l_0 = as<double>(par_0["loglikelihood"]);
  double Chisq = -2 * (l_0 - l_alt);

  if(as<int>(par_0["inflation"]) == 1){
    df = 2 * (p - 1);
  }
  else df = p - 1;

  // summarize results
  return List::create(_["par"] = par_alt, _["Chisq"] = Chisq, _["df"] = df, _["p.value"] = R::pchisq(Chisq, df, 0, 0));
}

// Multiple Testing
// [[Rcpp::export]]
List LRTnG(const arma::mat& X, const arma::mat& Y, const arma::vec& inflation, double tol = 1e-3, int maxIter = 5000, Nullable<List> initial = R_NilValue) {
  if (Y.n_cols == 1) {
    warning("It seems like you are testing a single variable. You may use `LRT1G`");
  }
  if (X.n_rows != Y.n_rows) {
    stop("The dimension of X and Y does not match!");
  }

  int k = Y.n_cols;

  arma::vec Chisq(k);
  arma::vec p_value(k);
  arma::vec df(k);

  for (int i = 0; i < k; ++i) {
    // keep running for multiple tests even if the error occurs
    try{
      List out = LRT1G(X, Y.col(i), tol, maxIter, inflation(i), initial);
      Chisq(i) = as<double>(out["Chisq"]);
      p_value(i) = as<double>(out["p.value"]);
      df(i) = as<double>(out["df"]);
    } catch(...) {
      warning("The regression model for one variable fails to converge, set its Chisq and p-value to be NA");
      Chisq(i) = NA_REAL;
      p_value(i) = NA_REAL;
      df(i) = NA_REAL;
    }
  }

  return List::create(_["Chisq"] = Chisq, _["df"] = df, _["p.value"] = p_value);
}
