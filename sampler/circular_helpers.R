#
# A small bag of helper functions required by circular sampler
# Pieter Jongsma 2015
#

require('circular')

Rcpp::sourceCpp('sampler/lib/rvmc.cpp')

#' Sample mean of von Mises distribution
SampleMu <- function(y, kappa) {
  if (length(y) == 0) {
    return(circular(runif(1, min=0, max=2*pi), template="clock24"))
  } else {
    estimation.kappa <- kappa*rho.circular(y)*length(y)
    return(circular(rvmc(1, mean(y), estimation.kappa), template="clock24"))
  }
}

#' Sample $\kappa$ according to Mardia scheme
#' Sampling can fail when observation count is small. In this case, NA is returned.
SampleKappa <- function(y, mu) {
  if (length(y) == 0) {
    return(NA)
  }
  
  R <- rho.circular(y)*length(y)
  etag <- -R * cos(mu - mean(y))
  eta <- length(y) # Uninformative
  
  res <- sampleKappa(etag, eta)
  if (res[1] < 0) {
    cat("Sampling kappa failed\n")
    return(NA)
  }
  
  return(res[1])
}

#' Makes sure all values fall within [mean-pi, mean+pi]
NormalizeMu <- function(mu) {
  mu[mu < mean(mu)-pi] <- mu[mu < mean(mu)-pi] + 2*pi
  mu[mu > mean(mu)+pi] <- mu[mu > mean(mu)+pi] - 2*pi
  mu
}

#' Converts radians to 24h clock range
HoursForRad <- function(mu) {
  1/(2*pi) * 24 * (as.numeric(mu) %% (2*pi))
}

#' String representatin for time
TimeForHours <- function(hours) {
  h <- as.integer(hours)
  m <- as.integer(hours*60 - h*60)
  paste0(str_pad(as.character(h), 2, "left", "0"), ":", str_pad(as.character(m), 2, "left", "0"))
}

#' String representatin for time
TimeForRad <- function(mu) {
  TimeForHours(HoursForRad(mu))
}
