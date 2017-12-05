#! /usr/local/bin/Rscript --vanilla

require(xtable)

# require(modeest) # Venter mode estimation
Rcpp::sourceCpp('sampler/lib/VenterMode.cpp')
Rcpp::sourceCpp('sampler/lib/VenterModeCirc.cpp')

source('configuration.R')
source('visualization.R')
source('sampler/circular_sampler.R')
source('helpers.R')

load("scenarios/three_components.RData")

pdf(file="figures/convergence.pdf", width=6, height=3)

indices <- 6006:6010
for (index in indices) {
  scenario <- scenarios[[index]]
  set.seed(scenario$seed)
  data <- circular(SimulateParams(scenario, n=scenario$n), template="clock24")
  load(paste0("results/", scenario$name, "-", index, ".RData"))

  output <- results
  ExpandOutput(output)

  log.likelihoods <- sapply(output, function(sample) {
    LogLikelihood(sample, data)
  })
  log.priorprobabilities <- sapply(output, function(sample) {
    LogPriorProbability(sample, data)
  })
  log.posterior <- log.likelihoods + log.likelihoods
  range <- 1:(burnin+iterations)

  par(mar=c(3, 5, 1, 1))
  if (index == indices[1]) {
    plot(range, log.posterior[range], type='l', col=custom.colors[[which(index == indices)]], ylim=c(-10^5, 0), xlab=NA, ylab="Log posterior probability")
  } else {
    lines(range, log.posterior[range], col=custom.colors[[which(index == indices)]])
  }
}

dev.off()



pdf(file="figures/convergence_zoomed.pdf", width=6, height=3)

indices <- 6006:6010
for (index in indices) {
  scenario <- scenarios[[index]]
  set.seed(scenario$seed)
  data <- circular(SimulateParams(scenario, n=scenario$n), template="clock24")
  load(paste0("results/", scenario$name, "-", index, ".RData"))

  output <- results
  ExpandOutput(output)

  log.likelihoods <- sapply(output, function(sample) {
    LogLikelihood(sample, data)
  })
  log.priorprobabilities <- sapply(output, function(sample) {
    LogPriorProbability(sample, data)
  })
  log.posterior <- log.likelihoods + log.likelihoods
  range <- (10^4):(1.5*10^4)

  par(mar=c(3, 5, 1, 1))
  if (index == indices[1]) {
    plot(range, log.posterior[range], type='l', col=custom.colors[[which(index == indices)]], ylim=c(-8000, -4000), xlab=NA, ylab="Log posterior probability")
  } else {
    lines(range, log.posterior[range], col=custom.colors[[which(index == indices)]])
  }
}

dev.off()
