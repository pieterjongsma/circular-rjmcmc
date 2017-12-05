
require('gtools')

ExpandOutput <- defmacro(output, expr={
  ks <- sapply(output, function(sample) {
    sample[['k']]
  })
  weights <- t(sapply(output, function(sample) {
    v <- rep(NA, max(ks))
    weight <- sample[['weight']]
    v[1:length(weight)] <- weight
    v
  }))
  mus <- t(sapply(output, function(sample) {
    v <- rep(NA, max(ks))
    mu <- sample[['mu']]
    # mu <- NormalizeMu(mu)
    v[1:length(mu)] <- mu
    v
  }))
  kappas <- t(sapply(output, function(sample) {
    v <- rep(NA, max(ks))
    kappa <- sample[['kappa']]
    v[1:length(kappa)] <- kappa
    v
  }))
})

SimulateParams <- function(params, n=100) {
  ExpandParameters(params)
  
  sapply(1:n, function(i) {
    j <- sample(1:k, 1, prob=weight)
    rvonmises(1, mu=mu[j], kappa=kappa[j])
  })
}

InitialParams <- function(data, k=1) {
  # Configuration
  mu <- circular(runif(k, min=0, max=2*pi), template="clock24")
  kappa <- runif(k, min=1, max=10)
  
  # Generated
  n <- length(data)
  weight <- rdirichlet(1, rep(n/k, k))
  z <- sample(1:k, n, replace=TRUE)
  
  list(mu=mu, kappa=kappa, k=k, z=z, weight=weight)
}
