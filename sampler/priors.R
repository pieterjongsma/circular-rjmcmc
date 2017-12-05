#' Calculate the prior probability of all parameters
LogPriorProbability <- function(p, data) {
  ExpandParameters(p)
  
  loglikelihood.weights <- sum(LogPriorWeight(weight))
  loglikelihood.mus <- sum(LogPriorMu(mu))
  
  loglikelihood.kappas <- sum(LogPriorKappa(kappa))
  
  loglikelihood.k <- LogPriorK(k, data)
  
  loglikelihood.alloc <- sum(LogPriorAllocation(p))

  # return(loglikelihood.weights + loglikelihood.mus + loglikelihood.kappas + loglikelihood.k)
  return(loglikelihood.alloc + loglikelihood.weights + loglikelihood.mus + loglikelihood.kappas + loglikelihood.k)
}
LogPriorAllocation <- function(params) {
  ExpandParameters(params)
  
  allocation.count <- sapply(1:k, function(j) sum(allocation == j))
  prob <- ddirichlet(weight, allocation.count+1)
  if (prob == 0) prob <- .Machine$double.eps
  log(prob)
  
  # allocation.count*log(weight)
}
LogPriorWeight <- function(weight) {
  log(ddirichlet(weight, rep(1, length(weight))))
  # dunif(weight, log=TRUE)
  # dbeta(weight, 3, 1, log=TRUE)
}
LogPriorMu <- function(mu) {
  dunif(mu %% 2*pi, min=0, max=2*pi, log=TRUE)
}
LogPriorKappa <- function(kappa) {
  0
}
LogPriorK <- function(k, data) {
  length(data) * dgeom(k, 0.05, log=TRUE)
}
PriorDrawMu <- function() {
  runif(1, min=0, max=2*pi)
}