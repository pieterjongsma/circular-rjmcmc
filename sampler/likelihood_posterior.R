#' Calculate the loglikelihood
#' This is the sum of the loglikelihoods of each observation for the component it is allocated to
LogLikelihood <- function(p, data) {
  ExpandParameters(p)
  
  likelihoods <- sapply(1:k, function(j) {
    data.subset <- data[allocation == j]
    sum(dvonmises(data.subset, mu=mu[j], kappa=kappa[j], log=TRUE))
  })
  
  return(sum(likelihoods))
}

LogPosteriorProbability <- function(params, data) {
  return(LogLikelihood(params, data) + LogPriorProbability(params, data))
}
