
weights_from_lik <- function(lik) {
  n <- length(lik)
  N <- 2^n
  ll <- log(lik)
  A <- as.matrix(expand.grid(lapply(1:n, \(i)(0:1)))) # 2^n x n
  stopifnot(is.matrix(A))
  stopifnot(nrow(A) == N)
  stopifnot(ncol(A) == n)
  kappa <- rowSums(A)
  stopifnot(length(kappa) == N)
  stopifnot(kappa %in% (0:n))
  u <- (A - 0.5) %*% ll
  stopifnot(length(u) == N)
  kvec <- 0:n
  # logw0 is log(w0) where w0 is the vector of unnormalized weights
  logw0 <-
    sapply(kvec, \(k){ logsumexp(u[kappa == k]) }) +
    lbeta(kvec + 1, n - kvec + 1)
  # logZ is the log of the normalizing constant
  logZ <- logsumexp(logw0)
  exp(logw0 - logZ)
}

logsumexp <- function(x) {
  xmax <- max(x)
  xmax + log(sum(exp(x - xmax)))
}

mean_from_weights <- function(w) {
  stopifnot(length(w) >= 1)
  n <- length(w) - 1
  k <- (0:n)
  alpha <- k + 1
  beta <- n - k + 1
  sum(w * alpha / (alpha + beta))
}

mean_from_lik <- function(lik) {
  mean_from_weights(weights_from_lik(lik))
}
