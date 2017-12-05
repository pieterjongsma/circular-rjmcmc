#' Update weights $\bm w$ according to allocation
AttemptMove.UpdateWeights <- function(params, data) {
  ExpandParameters(params)
  params.copy <- params
  
  allocation.count <- sapply(1:k, function(j) sum(allocation == j))
  new.weights <- rdirichlet(1, allocation.count+1)
  
  params.copy[['weight']] <- new.weights / sum(new.weights) # Assign normalized
  return(params.copy)
}

#' Estimate mu and kappa
AttemptMove.EstimateParameters <- function(params, data) {
  params <- AttemptMove.EstimateMu(params, data)
  params <- AttemptMove.EstimateKappa(params, data)
  params
}

#' Sample mu for each component using `SampleMu`
AttemptMove.EstimateMu <- function(params, data) {
  ExpandParameters(params)
  proposal <- params
  
  mu.new <- sapply(1:k, function(j) {
    y <- data[allocation == j]
    SampleMu(y, kappa[j])
  })
  
  mu.new <- mu.new[order((mu.new - 0.5*pi) %% (2*pi))]
  
  proposal[['mu']] <- circular(mu.new, template="clock24")
  
  return(proposal)
}

#' Sample kappa for each component using `SampleKappa` 
AttemptMove.EstimateKappa <- function(params, data) {
  ExpandParameters(params)
  proposal <- params
  
  default.kappa <- kappa
  
  kappa.new <- sapply(1:k, function(j) {
    y <- data[allocation == j]
    kap <- SampleKappa(y, mu[j])
    if (is.na(kap)) kap <- default.kappa[j] # Retain previous kappa when sampling fails
    kap
  })
  proposal[['kappa']] <- kappa.new
  
  return(proposal)
}

#' Update allocation according to the relative probability of an observation for each component
#' This method is quite inefficient
AttemptMove.UpdateAllocation <- function(params, data) {
  ExpandParameters(params)
  n <- length(data)
  params.copy <- params
  
  probs <- sapply(1:k, function(j) {
    weight[j] * dvonmises(data, mu=mu[j], kappa=kappa[j])
  })
  
  new.allocation <- row.apply(probs, function(prob) {
    prob[prob == 0 | is.na(prob)] <- .Machine$double.eps
    sample(1:k, 1, prob=prob/sum(prob))
  })
  
  params.copy[['z']] <- new.allocation
  return(params.copy)
}
