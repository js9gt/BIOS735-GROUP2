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

// [[Rcpp::export]]
arma::vec l1_theta(const arma::vec& z_k, const arma::vec& y_k, const arma::vec& mu_k, double thetak) {
  int n = y_k.n_elem;
  arma::vec dg(n);
  for (int i = 0; i < n; i++) {
    dg(i) = R::digamma(y_k(i) + 1 / thetak);
  }
  arma::vec score = (1 - z_k) / pow(thetak, 2) %
    (R::digamma(1 / thetak) - dg + log(1 + thetak * mu_k)) +
    (y_k - mu_k) % (1 - z_k) / thetak / (1 + thetak * mu_k);
  return score;
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

    arma::vec result(2 * p + 1);
    result.subvec(0, p - 1) = X.t() * V * ey;
    result.subvec(p, 2 * p - 1) = X.t() * ez;
    result(2 * p) = sum(l1_theta(z_k, y_k, mu_k, thetak)) * thetak;

    return result;
}

// Compute the approximate fisher information
// [[Rcpp::export]]
arma::vec l2_theta(const arma::vec& z_k, const arma::vec& y_k, const arma::vec& mu_k, double thetak) {
  int n = y_k.n_elem;
  arma::vec dg(n);
  arma::vec tg(n);
  for (int i = 0; i < n; i++) {
    dg(i) = R::digamma(y_k(i) + 1 / thetak);
    tg(i) = R::trigamma(y_k(i) + 1 / thetak);
  }
  arma::vec l2 = (1 - z_k) / pow(thetak, 4) % (
    2 * thetak * R::digamma(1 / thetak) - 2 * thetak * dg +
    2 * thetak * (1 + thetak * mu_k) + R::trigamma(1 / thetak) - tg -
    2 * pow(thetak, 2) * mu_k / (1 + thetak) -
    2 * pow(thetak, 2) * mu_k / (1 + thetak * mu_k)
  );
  return l2;
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

  arma::mat I_mat = arma::zeros<arma::mat>(2 * p + 1, 2 * p + 1);
  I_mat.submat(0, 0, p - 1, p - 1) = X.t() * Py * X;
  I_mat.submat(p, p, 2 * p - 1, 2 * p - 1) = X.t() * Pz * X;
  I_mat(2 * p, 2 * p) = sum(l2_theta(z_k, y_k, mu_k, thetak)) * pow(thetak, 2) + sum(l1_theta(z_k, y_k, mu_k, thetak)) * thetak;

  return I_mat;
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
  long double thetak, loglikelihood;

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

  // determin whether inflated model should be fitted
  if(min(y_k) != 0) {
    inflation = false;
  }

  int iter = 0;
  double eps = INFINITY;
  arma::mat I_mat, identity;

  // zero inflated model
  if(inflation){
    identity = arma::eye(2 * p + 1, 2 * p + 1);
    while (iter < maxIter && eps > tol) {
      // E step
      arma::vec pi_k = c_pi(X, gammak);
      arma::vec mu_k = c_mu(X, betak);
      arma::vec z_k = Qik(pi_k, mu_k, y_k, thetak);

      // M step using one step Newton-Raphson method
      arma::vec s = score(X, y_k, z_k, mu_k, pi_k, thetak);
      I_mat = In(X, y_k, z_k, mu_k, pi_k, gammak, thetak);

      arma::vec par_old = join_vert(join_vert(betak, gammak), arma::vec({log(thetak)}));
      arma::vec par_new = par_old + arma::solve(I_mat, s, arma::solve_opts::fast + arma::solve_opts::allow_ugly);
      eps = sqrt(sum(square(par_old - par_new)));

      betak = par_new.subvec(0, p - 1);
      gammak = par_new.subvec(p, 2 * p - 1);
      thetak = exp(par_new(2 * p));

      iter = iter + 1;
    }
    // compute the loglikelihood
    arma::vec pi_k = c_pi(X, gammak);
    arma::vec mu_k = c_mu(X, betak);
    arma::vec p_0 = 1 / (1 + thetak * mu_k);
    arma::vec logprob = log(dznbinom(y_k, 1 / thetak, p_0, pi_k));
    // remove NA and Infs
    loglikelihood = accu(logprob.elem(find_finite(logprob)));
  }
  else{
    // this is the given information for the non-inflated model
    identity = arma::eye(p + 1, p + 1);
    arma::vec pi_k = arma::zeros(n);
    arma::vec z_k = arma::zeros(n);
    arma::uvec index = arma::join_cols(arma::regspace<arma::uvec>(0, p - 1), arma::uvec({static_cast<unsigned int>(2 * p)}));

    // fit a negative binomial regression
    while (iter < maxIter && eps > tol) {
      arma::vec mu_k = c_mu(X, betak);
      // one step Newton-Raphson method
      arma::vec s = score(X, y_k, z_k, mu_k, pi_k, thetak)(index);
      I_mat = In(X, y_k, z_k, mu_k, pi_k, gammak, thetak)(index, index);

      arma::vec par_old = join_vert(betak, arma::vec({log(thetak)}));
      arma::vec par_new = par_old + arma::solve(I_mat, s, arma::solve_opts::fast + arma::solve_opts::allow_ugly);
      eps = sqrt(sum(square(par_old - par_new)));

      betak = par_new.subvec(0, p - 1);
      thetak = exp(par_new(p));
      iter = iter + 1;
    }
    // compute the loglikelihood
    arma::vec mu_k = c_mu(X, betak);
    arma::vec p_0 = 1 / (1 + thetak * mu_k);
    arma::vec logprob = log(dznbinom(y_k, 1 / thetak, p_0, pi_k));
    // remove NA and Infs
    loglikelihood = accu(logprob.elem(find_finite(logprob)));
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
    _["loglikelihood"] = loglikelihood
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
