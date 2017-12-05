#! /usr/local/bin/Rscript --vanilla

source('configuration.R')
source('sampler/circular_sampler.R')

EvaluateScenarios <- function(input.file, output.file, start.index=1) {
  load(input.file) # Creates variable `scenarios`
  
  for (index in start.index:length(scenarios)) {
    scenario <- scenarios[[index]]
  
    results <- EvaluateScenario(scenario)
  
    output <- paste0(output.file, "-", index, ".RData")
    save(results, file=output)
  }
}

EvaluateScenario <- function(scenario) {
  set.seed(scenario$seed)
  data <- circular(SimulateParams(scenario, n=scenario$n), template="clock24")
  
  start <- InitialParams(data)
  results <- RunSampler(data, start, burnin+iterations)
  
  return(results)
}


args <- commandArgs(TRUE)
input.file <- args[1]
output.file <- args[2]
start.index <- as.integer(args[3])
if (is.na(start.index)) { start.index <- 1 }

load('plays.RData')
EvaluateScenarios(input.file, output.file, start.index)
