#! /usr/local/bin/Rscript --vanilla

#
# Script generating csv file with instructions for simulations
# Pieter Jongsma 2015
#

source('configuration.R')

GenerateSeed <- function(n) {
  sample(1:10^9, n)
}

GenerateScenarios <- function(definition, replications=1) {
  ns <- c(5 * 10^4)
  
  scenarios <- lapply(ns, function(n) {
    modifyList(definition, list(n=n))
  })
  scenarios.rep <- scenarios[rep(1:length(scenarios), each=replications)]
  
  seeds <- GenerateSeed(length(scenarios.rep))
  scenarios.rep <- lapply(1:length(scenarios.rep), function(i) {
    modifyList(scenarios.rep[[i]], list(seed=seeds[i]))
  })
  
  scenarios.rep
}

for (definition in scenario.definitions) {
  scenarios <- GenerateScenarios(definition, 1000)
  file <- paste0('scenarios_big/', definition$name, '.RData')
  save(scenarios, file=file)
}
cat("Generated scenarios.\n")
