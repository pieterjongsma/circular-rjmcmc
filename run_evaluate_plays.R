
require('circular')
require('Rcpp')
require('stringr')

source('configuration.R')
source('visualization.R')
source('sampler/circular_sampler.R')

Rcpp::sourceCpp('sampler/lib/VenterMode.cpp')

set.seed(18081992)

load('plays.RData')
plays.counts <- sapply(plays, length)

burnin <- 10^4
iterations <- 10^4


data <- circular(unlist(plays), template="clock24")
start <- InitialParams(data)
results <- RunSampler(data, start, burnin+iterations)

output <- results[(burnin+1):(burnin+iterations)]

SummarizeOutput <- function(output) {
  ExpandOutput(output)

  mus <- t(sapply(output, function(sample) {
    v <- rep(NA, max(ks))
    mu <- sample$mu
    mu <- NormalizeMu(mu)
    v[1:length(mu)] <- mu
    v
  }))
  kappas <- t(sapply(output, function(sample) {
    v <- rep(NA, max(ks))
    kappa <- sample$kappa
    v[1:length(kappa)] <- kappa
    v
  }))
  weights <- t(sapply(output, function(sample) {
    v <- rep(NA, max(ks)+1)
    weight <- sample$weight
    v[1:length(weight)] <- weight
    v
  }))
  if (max(ks) > 1) {
    permutation <- t(row.apply(mus, function(mu) {
      order(mu)
    }))
    for (i in 1:nrow(permutation)) {
      mus[i, ] <- mus[i, permutation[i, ]]
      kappas[i, ] <- kappas[i, permutation[i, ]]
      weights[i, ] <- weights[i, permutation[i, ]]
    }
  } else {
    mus <- t(mus)
    kappas <- t(kappas)
  }

  k.mode <- unique(ks)[which.max(table(ks))]
  mu.mean <- circular(col.apply(as.matrix(mus[ks==k.mode, 1:k.mode]), mean), template="clock24")
  HoursForRad(mu.mean)

  kappa.mean <- col.apply(as.matrix(kappas[ks==k.mode, 1:k.mode]), function(x) {
    hmode(x, 0.10)
  })
  kappa.mean

  weight.mean <- col.apply(as.matrix(weights[ks==k.mode, 1:(k.mode+1)]), mean)
  weight.mean

  frame <- 1.8
  sample <- output[[10^4]]
  par(mfrow=c(1, 2), mar=c(0, 0, 0, 0))
  plot(sample$mu, xlim=c(-frame, frame), ylim=c(-frame, frame))
  CurvesForParameters(sample, add=TRUE)
  plot(sample$mu, xlim=c(-frame, frame), ylim=c(-frame, frame))
  CurveForParameters(sample, add=TRUE)

  scenario <- list(
    mu=mu.mean,
    kappa=kappa.mean,
    weight=weight.mean,
    allocation=sample(1:k.mode, length(data), replace=TRUE),
    k=k.mode
  )
  
  print(table(ks)/length(ks))
  for (j in 1:(scenario$k+1)) {
    cat("Weight :", scenario$weight[j])
    cat("; Mu :", TimeForRad(scenario$mu[j]))
    cat("; Kappa :", scenario$kappa[j])
    cat("; LB :", TimeForRad(qvonmises(0.025, mu=scenario$mu[j], kappa=scenario$kappa[j])))
    cat("; RB :", TimeForRad(qvonmises(0.975, mu=scenario$mu[j], kappa=scenario$kappa[j])))
    cat("\n")
  }
}

SummarizeOutput(output)

#
# Separate genres
#

source('sampler/circular_sampler_hacked.R')

load('plays_joined_scenario.RData')

set.seed(18081992)

plays.index <- 2
data <- plays[[plays.index]]
start <- InitialParams(data)
start$weight <- c(start$weight, 0.0)
results <- RunSampler(data, start, burnin+iterations, scenario)

output <- results[(burnin+1):(burnin+iterations)]
SummarizeOutput(output)
