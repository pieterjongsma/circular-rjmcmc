
prob <- function(i, j, theta) {
  x <- data[i]
  probs <- sapply(1:k, function(s) {
    theta$weight[s] * dvonmises(x, mu=theta$mu[s], kappa=theta$kappa[s])
  })
  prob <- probs[j] / sum(probs)
  return(prob)
}
apply.permutation <- function(theta, permutation) {
  list(
    mu=theta$mu[permutation],
    kappa=theta$kappa[permutation],
    weight=theta$weight[permutation]
    # FIXME Include z / allocation
    )
}
apply.permutation.list <- function(theta, permutation) {
  lapply(1:(length(theta)-1), function(t) { # FIXME Why length(theta)-1 ?
    th <- theta[[t]]
    apply.permutation(th, permutation[t, ])
  })
}

LS.UpdateAssignmentProb <- function(samples, permutation, k) { # Step 1
  n <- length(samples[[1]]$z)
  N <- length(samples)
  
  vector <- sapply(1:n, function(i) {
    sapply(1:k, function(j) {
      # Expression (17) Stephens
      som <- sum(sapply(1:N, function(t) {
        theta <- samples[[t]]
        theta <- apply.permutation(theta, permutation[t, ])
        prob(i, j, theta)
      }))
      som/N
    })
  })
  
  q <- t(vector)
  
  q
}

LS.UpdatePermutation <- function(samples, q, k) { # Step 2
  n <- length(samples[[1]]$z)
  N <- length(samples)
  selected.permutation = matrix(rep(NA, k*N), ncol=k)
  
  for (t in 1:N) {
    lowest.som <- -1
    for (permutation in permn(1:k)) {
      # Expression (16) Stephens
      som <- 0
      for (i in 1:n) {
        for (j in 1:k) {
          theta <- samples[[t]]
          theta <- apply.permutation(theta, permutation)
          p <- prob(i, j, theta)
          val <- p * log(p/q[i, j])
          som <- som + val
        }
      }
      if (lowest.som < 0 || som < lowest.som) {
        selected.permutation[t, ] <- permutation
        lowest.som <- som
      }
    }
  }
  
  selected.permutation
}

PostProcessLabels <- function(samples) {
  k <- samples[[1]]$k
  N <- length(samples)
  n <- length(samples[[1]]$z)

  permutation <- matrix(rep(1:k, each=N), ncol=k)
  permutation.prev <- NA
  for (repetition in 1:100) {
    cat("Repetition", repetition, "\n")
    q <- LS.UpdateAssignmentProb(samples, permutation, k)
    cat("Step 1 done\n")
    permutation <- LS.UpdatePermutation(samples, q, k)
    cat("Step 2 done\n")
  
    if (repetition > 1) {
      diff <- sum(permutation != permutation.prev)
      cat("Diff:", diff, "\n")
      if (diff == 0) {
        break
      }
    }
    permutation.prev <- permutation
  }
  
  return(apply.permutation.list(samples), permutation)
}
