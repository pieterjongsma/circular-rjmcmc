#
# Sampler implementation for fitting a mixture of von Mises distributions with a unknown component count
# Pieter Jongsma 2015
#

require('circular')
require('gtools') # defmacro, dirichlet functions

source('sampler/circular_helpers.R')
Rcpp::sourceCpp('sampler/lib/sampleKappa.cpp')

# Nice-ify
row.apply <- function(X, ...) { apply(X, 1, ...) }
col.apply <- function(X, ...) { apply(X, 2, ...) }

#' Function that puts parameters in params list into their own variables for easy access
ExpandParameters <- defmacro(p, expr={ # Shorthand for assigning all parameters in list object to named variables for easy referencing
  k <- p[['k']]; mu <- circular(p[['mu']], template="clock24"); kappa <- p[['kappa']]; allocation <- p[['z']]; weight <- p[['weight']]
})

# These files require the `ExpandParameters` method
source('sampler/priors.R')
source('sampler/likelihood_posterior.R')
source('sampler/fixed_dimension_moves.R')
source('sampler/dimensionality_changing_moves.R')

#' Run the RJMCMC sampler on the data
#' Returns the MCMC chain of all parameters as a list of lists
RunSampler <- function(data, startvalue, iterations, change.dimensions=TRUE) {
  chain <- list()
  
  chain[[1]] <- startvalue
  for (i in 1:iterations) {
    params <- chain[[i]]
    if (i %% 50 == 0) {
      cat("Iteration", i, "( k =", params[['k']], ")\n")
    }
    
    params <- AttemptMove.UpdateWeights(params, data)
    params <- AttemptMove.EstimateParameters(params, data)
    if (change.dimensions && i %% 5 == 0) {
      params <- AttemptMove.SplitCombine(params, data)
      params <- AttemptMove.BirthDeath(params, data)
    }
    params <- AttemptMove.UpdateAllocation(params, data)
    
    chain[[i+1]] <- params
  }
  
  return(chain)
}

log.info <- function(...) {
  if (LOG) { # Ability to switch to boolean for verbose / non-verbose
    cat(...)
  }
}
