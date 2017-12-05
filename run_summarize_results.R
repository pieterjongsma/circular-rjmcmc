#! /usr/local/bin/Rscript --vanilla

require(xtable)

# require(modeest) # Venter mode estimation
Rcpp::sourceCpp('sampler/lib/VenterMode.cpp')
Rcpp::sourceCpp('sampler/lib/VenterModeCirc.cpp')

source('configuration.R')
burnin <- 10^4 + 4*10^3
iterations <- 10^3
source('sampler/circular_helpers.R')
source('helpers.R')

row.apply <- function(X, ...) { apply(X, 1, ...) }
col.apply <- function(X, ...) { apply(X, 2, ...) }

IsInIntervalCirc <- function(true.value, interval) {
  if (interval[1] < interval[2] && (true.value > interval[1] && true.value < interval[2])) {
    1
  } else if (interval[1] > interval[2] && (true.value > interval[1] || true.value < interval[2])) {
    1
  } else {
    0
  }
}



scenario.indices <- 1:length(scenario.definitions)
scenario.indices <- 1:5
results.dir <- "results/"
summaries <- list()

for (scenario.i in scenario.indices) {
  cat("Summarizing scenario group", scenario.i, "\n")
  
  definition <- scenario.definitions[[scenario.i]]
  scenario.file <- paste0('scenarios/', definition$name, '.RData')
  load(scenario.file) # Creates `scenarios`

  ns <- unique(sapply(scenarios, function(s) {
    s$n
  }))

  ks.counts <- list()
  mu.estimates <- matrix(nrow=length(scenarios), ncol=definition$k)
  kappa.estimates <- mu.estimates
  weight.estimates <- mu.estimates

  for (i in 1:length(scenarios)) {
  # for (i in 1:5) {
    cat("Summarizing scenario", i, "\n")
    
    results.file <- paste0(results.dir, definition$name, "-", i, ".RData")
    if (file.exists(results.file)) {
      load(results.file) # Creates `results`
      sample.range <- (burnin+1):(burnin+iterations)
      samples <- results[sample.range]
      ExpandOutput(samples)
    
      k <- definition$k
    
      ks.u <- 1:max(ks)
      ks.count <- sapply(ks.u, function(j) { sum(ks == j) })
      ks.counts[[i]] <- ks.count
    
    
      subset <- (ks == definition$k)
      samples.subset <- samples[subset]
      ExpandOutput(samples.subset)
      
      if (ncol(mus) < k) { next }
    
      permutation <- t(row.apply(as.matrix(mus[, 1:k]), function(mu) {
        order((mu - pi/2) %% (2*pi))
      }))
      mus <- t(sapply(1:nrow(mus), function(i) {
        r <- mus[i, permutation[i, ]]
        r %% (2*pi)
      }))
      kappas <- t(sapply(1:nrow(kappas), function(i) {
        r <- kappas[i, permutation[i, ]]
      }))
      weights <- t(sapply(1:nrow(weights), function(i) {
        r <- weights[i, permutation[i, ]]
      }))
      
      mu.estimates[i, ] <- sapply(1:definition$k, function(j) {
        estimate <- hmodecirc(mus[, j], 0.1)
        return(estimate)
      })
      kappa.estimates[i, ] <- sapply(1:definition$k, function(j) {
        estimate <- hmode(kappas[, j], 0.1)
        return(estimate)
      })
      weight.estimates[i, ] <- sapply(1:definition$k, function(j) {
        estimate <- hmode(weights[, j], 0.1)
        return(estimate)
      })
    }
  }
  
  summaries[[scenario.i]] <- list(ks.counts=ks.counts, mu.estimates=mu.estimates, kappa.estimates=kappa.estimates, weight.estimates=weight.estimates)
}


for (scenario.i in scenario.indices) {
  definition <- scenario.definitions[[scenario.i]]
  summary <- summaries[[scenario.i]]
  
  ks.counts <- summary$ks.counts
  mu.estimates <- summary$mu.estimates
  
  k.modes <- lapply(ks.counts, function(counts) {
    which.max(counts)
  })
  
  freqs <- sapply(ks.counts, function(counts) {
    cutoff <- 5
    if (length(counts) == 0) { return(rep(NA, cutoff)) }
    freq <- counts / sum(counts)
    freq <- c(freq[1:(cutoff-1)], sum(freq[cutoff:max(cutoff, length(freq))]))
    freq[is.na(freq)] <- 0
    return(freq)
  })

  frame <- data.frame()
  for (n in ns) {
    indices <- which(sapply(scenarios, function(s) { s$n == n }))
    
    indices <- indices[indices < length(ks.counts)]
    
    cat("Scenario", definition$name, ", n =", n, ":\n")
    
    modes <- unlist(k.modes[indices])
    correct.mode <- mean(modes == definition$k)
    cat("    k mode correct (0-1) :", correct.mode, "\n")
    
    freqs.mean <- row.apply(as.matrix(freqs[, indices]), function(x) mean(x, na.rm=TRUE))
    cat("    k freqs :", freqs.mean, "\n")
    
    for (j in 1:definition$k) {
      estimate <- mean(circular(mu.estimates[indices, j], template="clock24"), na.rm=TRUE)
      estimate <- estimate[!is.na(estimate)]
      estimate[estimate > pi] <- estimate[estimate > pi] - 2*pi
      cat("    mu", j, ":", estimate, "\n")
      
      estimate <- mean(kappa.estimates[indices, j], na.rm=TRUE)
      estimate <- estimate[!is.na(estimate)]
      cat("    kappa", j, ":", estimate, "\n")
      
      estimate <- mean(weight.estimates[indices, j], na.rm=TRUE)
      estimate <- estimate[!is.na(estimate)]
      cat("    weight", j, ":", estimate, "\n")
    }
    
    cat("    Count :", length(modes), "\n")
    
    frame <- rbind(frame, data.frame(n=n, k.correct=correct.mode, freqs=t(freqs.mean), replications=length(modes)))
  }
  
  print(xtable(frame, digits=c(0, 0, 2, 2, 2, 2, 2, 2, 0)))
}
