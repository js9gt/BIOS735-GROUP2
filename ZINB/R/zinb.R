# E step compute the Q function given parameters in the previous step
c.pi <- function(X, gammak){
  exp(X %*% gammak) / (1 + exp(X %*% gammak))
}

c.mu <- function(X, betak){
  exp(X %*% betak)
}

Qik <- function(piik, muik, yik, thetak){
  indicator = as.numeric(yik == 0)
  return(piik * indicator /
           (piik + (1 - piik) / (1 + thetak * muik) ^ (1/thetak)))
}

# M step helper functions
# Score with respect to theta
l1_theta <- function(zik, yik, muik, thetak){
  score = (1 - zik) / thetak ^ 2 *
    (digamma(1/thetak) - digamma(yik+1/thetak) + log(1+thetak*muik)) +
    (yik - muik) * (1 - zik) / thetak / (1 + thetak * muik)
  return(score)
}

# Score function: Compute the score of the function
score <- function(X, y_k, z_k, mu_k, pi_k, thetak){
  n = dim(X)[1]

  V = matrix(0, n, n)
  diag(V) = 1 / (1 + thetak * mu_k)
  ey = (1 - z_k) * (y_k - mu_k)
  ez = z_k - pi_k
  return(
    c(
      t(X) %*% V %*% ey,
      t(X) %*% ez,
      sum(l1_theta(z_k, y_k, mu_k, thetak))
    )
  )
}

# Compute the approximate fisher information
l2_theta <- function(zik, yik, muik, thetak){
  l2 = (1 - zik) / thetak ^ 4 * (
    2 * thetak * digamma(1/thetak) - 2 * thetak * digamma(yik+1/thetak) +
      2 * thetak * (1 + thetak * muik) + trigamma(1/thetak) -
      trigamma(yik + 1 / thetak) -
      2 * thetak^2 * muik / (1 + thetak * muik)
  )

  return(l2)
}

In <- function(X, y_k, z_k, mu_k, pi_k, gammak, thetak){
  n = dim(X)[1]
  p = dim(X)[2]

  Py = Pz = matrix(0, n, n)
  diag(Py) = (1 - z_k) * mu_k / (1 + thetak * mu_k)
  diag(Pz) = pi_k / (1 + exp(X %*% gammak))

  I_mat = matrix(0, 2*p+1, 2*p+1)
  I_mat[1:p, 1:p] = t(X) %*% Py %*% X
  I_mat[(p+1):(2*p), (p+1):(2*p)] = t(X) %*% Pz %*% X
  I_mat[2*p+1, 2*p+1] = sum(l2_theta(z_k, y_k, mu_k, thetak))

  return(I_mat)
}

# EM algorithm
par.EM <- function(X, y_k, tol=1e-4, maxIter=1000, initial=NULL){
  n = dim(X)[1]
  p = dim(X)[2]

  # you can specify the initial points of beta, gamma and theta
  if(is.null(initial)){
    betak = gammak = rnorm(2)
    thetak = runif(1, min=0, max=1)
  }
  else if(length(initial) == 3){
    betak = initial[[1]]
    gammak = initial[[2]]
    thetak = initial[[3]]
  }
  iter = 1
  eps = Inf
  while(iter <= 1000 & eps > tol){
    # E step
    pi_k = c.pi(X, gammak)
    mu_k = c.mu(X, betak)
    z_k = Qik(pi_k, mu_k, y_k, thetak)
    p_0 = 1 / (1 + thetak * mu_k)
    print(log(prod(dnbinom(x=y_k, size=1/thetak, prob=p_0))))

    # M step using one step NR
    s = score(X, y_k, z_k, mu_k, pi_k, thetak)
    I_mat = In(X, y_k, z_k, mu_k, pi_k, gammak, thetak)
    print(I_mat)

    par_old = c(betak, gammak, thetak)
    par_new = par_old + solve(I_mat) %*% s
    eps = sqrt(sum((par_old - par_new)^2) / (sum(par_old^2)))

    betak = par_new[1:p]
    gammak = par_new[(p+1):(2*p)]
    thetak = par_new[2*p+1]
    print(par_new)

    iter = iter + 1
  }

  return(list(beta = betak, gamma = gammak, theta = thetak))
}

# Test Significance for single species
LRT1G <- function(X, y_k){
  n = dim(X)[1]
  p = dim(X)[2]

  # alternative hypothesis
  par_alt = par.EM(X, y_k)
  mu_k = c.mu(X, par_alt$beta)
  p = 1 / (1 + par_alt$theta * mu_k)
  L_alt = prod(dnbinom(x=y_k, size=1/par_alt$theta, prob=p))

  # null hypothesis
  X_0 = matrix(rep(0, n), n, 1)
  par_0 = par.EM(X_0, y_k)
  mu_k = c.mu(X, par_0$beta)
  p = 1 / (1 + par_0$theta * mu_k)
  L_0 = prod(dnbinom(x=y_k, size=1/par_0$theta, prob=p))
  Chisq = -2 * (log(L_0) - log(L_alt))

  # summarize results
  list(
    par = par_alt,
    loglikelihood = log(L_alt),
    Chisq = Chisq,
    df = 2 * (p - 1),
    p.value = 1 - pchisq(Chisq, 2 * (p - 1))
  )
}
