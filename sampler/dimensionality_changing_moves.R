#' Returns either the original parameters or the new proposal based on a Metropolis Hastings ratio
#' Has some additional error handling for when the ratio is NA or NaN
EvalMHRatio <- function(mh.ratio, params, proposal) {
  log.info("MH ratio", mh.ratio, "\n")
  if (!is.na(mh.ratio) && !is.nan(mh.ratio) && runif(1) < mh.ratio) {
    return(proposal)
  } else {
    return(params)
  }
}

#' Generate helper vector $\bm u$ that is used to create a split proposal
GenerateSplitCombineU <- function() {
  u1 <- 0.5 * runif(1)
  u2 <- runif(1, min=0, max=2*pi)
  u3 <- rbeta(1, 2, 1)
  c(u1, u2, u3)
}

#' Calculate the logprob of the helper vector $\bm u$
LogProbabilitySplitCombineU <- function(u) {
  prob <- c(
    dunif(u[1], log=TRUE),
    dunif(u[2], min=0, max=2*pi, log=TRUE),
    dbeta(u[3], 2, 1, log=TRUE)
    )
  sum(prob)
}

#' Perform a Metropolis Hastings split move 
#' 1. Generate helper vector $\bm u$
#' 2. Calculate the split result from this helper vector
AttemptMove.Split <- function(params, data) {
  ExpandParameters(params)
  
  log.info("Split ( k =", k, ")\n") # Will be appended to with MH ratio
  
  log.move.prob.ratio <- log(SplitCombine.MoveRatio(params, data))
  
  u <- GenerateSplitCombineU()
  log.u.prob <- LogProbabilitySplitCombineU(u)
  
  j <- sample(1:k, 1)
  weight.j <- weight[j]
  mu.j <- mu[j]
  kappa.j <- kappa[j]
  
  R.j <- A1(kappa.j)
  
  gamma <- u[2]
  rho <- (1-R.j) * u[3]
  
  weight.1 <- weight.j * u[1]
  weight.2 <- weight.j - weight.1
  
  a.1 <- R.j - cos(gamma) * rho
  b.1 <- sin(gamma) * rho
  R.1 <- sqrt(a.1^2 + b.1^2)
  mu.1 <- circular(mu.j - atan2(b.1, a.1), template="clock24")
  
  a.2 <- R.j - cos(gamma+pi) * rho * weight.1/weight.2
  b.2 <- sin(gamma+pi) * rho * weight.1/weight.2
  R.2 <- sqrt(a.2^2 + b.2^2)
  mu.2 <- circular(mu.j - atan2(b.2, a.2), template="clock24")
  
  kappa.1 <- A1inv(R.1); kappa.2 <- A1inv(R.2)
  
  weight[j] <- weight.1
  weight <- append(weight, weight.2, after=j)
  mu[j] <- mu.1
  mu <- append(mu, mu.2, after=j)
  kappa[j] <- kappa.1
  kappa <- append(kappa, kappa.2, after=j)
  
  log.info("Splitting weight", weight.j, "into", weight.1, "and", weight.2, "\n")
  log.info("Splitting mu", mu.j, "into", mu.1, "and", mu.2, "\n")
  log.info("Splitting kappa", kappa.j, "into", kappa.1, "and", kappa.2, "\n")
  
  proposal <- modifyList(params, list(weight=weight, mu=mu, kappa=kappa, k=k+1))
  
  # Update allocation for split
  idx <- c(j, j+1)
  const <- 1/(2*pi*besselI(kappa[idx], 0)) # These values take long to calculate, and are constant for every component, so cache
  const[is.na(const)] <- 0
  new.allocation <- sapply(data[allocation == j], function(y) {
    prob <- as.double(weight[idx] * const * exp(kappa[idx] * cos(y-mu[idx])))
    # prob <- c(weight.1, weight.2)
    prob <- prob/sum(prob)
    prob[prob == 0 | is.na(prob)] <- .Machine$double.eps
    sample(idx, 1, prob=prob)
  })
  if (j < k) {
    allocation[allocation > j] <- allocation[allocation > j] + 1
  }
  if (length(new.allocation) > 0) {
    allocation[allocation == j] <- new.allocation
  }
  
  proposal <- modifyList(proposal, list(z=allocation))
  
  comb.r <- R.j
  comb.weight <- weight.j
  jacobian <- abs ( ((comb.r - 1)^2*comb.r*comb.weight*(1-2*u[1])^2*u[3]) / (
      (u[1]-1)^2 *
      sqrt(2*(comb.r-1) * comb.r * cos(u[2]) * u[3] + u[3]^2 - 2*comb.r*u[3]^2 + comb.r^2*(1+u[3]^2)) *
      sqrt( 1/(u[1]-1)^2 * (-2*(comb.r-1)*comb.r*cos(u[2])*(u[1]-1)*u[1]*u[3] + u[1]^2*u[3]^2 - 2*comb.r*u[1]^2*u[3]^2 + comb.r^2*(1-2*u[1]+u[1]^2*(1+u[3]^2))) )
      ) )
  
  log.info("Split jacobian", jacobian, "\n")

  log.proposal <- LogPosteriorProbability(proposal, data)
  log.params <- LogPosteriorProbability(params, data)
  A <- exp(log.proposal - log.params + log.move.prob.ratio - log.u.prob + log(jacobian))
  log.info("Proposal prob:", exp(log.proposal), "(log:", log.proposal, ")", "\nParams prob:", exp(log.params), "(log:", log.params, ")\n")
  mh.ratio <- min(1, A)
  
  EvalMHRatio(mh.ratio, params, proposal)
}
AttemptMove.Combine <- function(params, data) {
  ExpandParameters(params)
  
  log.info("Combine ( k =", k, ")\n")
  
  if (k <= 1) {
    return(params)
  }
  
  log.move.prob.ratio <- log(SplitCombine.MoveRatio(params, data))
  
  j <- sample(1:k, 1)
  # jj <- sample((1:k)[-j], 1)
  jj <- order(abs(mu-mu[j]))[2] # Find the closest mu
  idx <- c(j, jj)
  
  comb.weight <- weight[j] + weight[jj]
  u1 <- weight[j] / comb.weight
  
  r.j <- A1(kappa[j]); r.jj <- A1(kappa[jj])
  a.j <- cos(mu[j]) * r.j
  b.j <- sin(mu[j]) * r.j
  a.jj <- cos(mu[jj]) * r.jj
  b.jj <- sin(mu[jj]) * r.jj
  
  rho.a <- (a.jj - a.j) * (weight[jj] / comb.weight)
  rho.b <- (b.jj - b.j) * (weight[jj] / comb.weight)
  rho <- sqrt(rho.a^2 + rho.b^2)
  comb.a <- a.j + rho.a
  comb.b <- b.j + rho.b
  comb.mu <- circular(atan2(comb.b, comb.a), template="clock24")
  comb.r <- sqrt(comb.a^2 + comb.b^2)
  
  b <- r.j - comb.r * cos(comb.mu - mu[jj])
  s <- sqrt((comb.a-a.j)^2 + (comb.b-b.j)^2)
  gamma <- runif(1, min=0, max=2*pi) # FIXME Actually derive gamma. Right now though, it does not influence probability, because of uniform prior for mu
  
  u2 <- gamma %% (2*pi)
  u3 <- rho / (1 - comb.r)
  
  comb.kappa <- A1inv(comb.r)
  
  new.allocation <- allocation
  new.allocation[new.allocation == max(idx)] <- min(idx)
  new.allocation[new.allocation > max(idx)] <- new.allocation[new.allocation > max(idx)]-1
  
  log.info("Combining weight", weight[jj], "and", weight[j], "into", comb.weight, "\n")
  log.info("Combining mu", mu[jj],         "and", mu[j],     "into", comb.mu, "\n")
  log.info("Combining kappa", kappa[jj],   "and", kappa[j],  "into", comb.kappa, "\n")
  
  weight[j] <- comb.weight
  weight <- weight[-jj]
  mu[j] <- comb.mu
  mu <- mu[-jj]
  kappa[j] <- comb.kappa
  kappa <- kappa[-jj]
  
  proposal <- modifyList(params, list(weight=weight, mu=mu, kappa=kappa, z=new.allocation, k=k-1))
  
  u <- c(u1, u2, u3)
  log.u.prob <- LogProbabilitySplitCombineU(u)
  
  jacobian <- abs ( ((comb.r - 1)^2*comb.r*comb.weight*(1-2*u[1])^2*u[3]) / (
      (u[1]-1)^2 *
      sqrt(2*(comb.r-1) * comb.r * cos(u[2]) * u[3] + u[3]^2 - 2*comb.r*u[3]^2 + comb.r^2*(1+u[3]^2)) *
      sqrt( 1/(u[1]-1)^2 * (-2*(comb.r-1)*comb.r*cos(u[2])*(u[1]-1)*u[1]*u[3] + u[1]^2*u[3]^2 - 2*comb.r*u[1]^2*u[3]^2 + comb.r^2*(1-2*u[1]+u[1]^2*(1+u[3]^2))) )
      ) )
  
  log.proposal <- LogPosteriorProbability(proposal, data)
  log.params <- LogPosteriorProbability(params, data)
  A <- exp(log.proposal - log.params - log.move.prob.ratio + log.u.prob - log(jacobian)) # Note that this A is the inverse of split move A
  log.info("Proposal prob:", exp(log.proposal), "(log:", log.proposal, ")", "\nParams prob:", exp(log.params), "(log:", log.params, ")\n")
  mh.ratio <- min(1, A)
  EvalMHRatio(mh.ratio, params, proposal)
}
SplitCombine.MoveRatio <- function(params, data) {
  ExpandParameters(params)
  
  bk <- 0.5 # Prob of combine # FIXME Shouldn't this be 0 if k=1
  dk1 <- (1-bk) # Prob of split
  return(dk1 / bk)
}
AttemptMove.SplitCombine <- function(params, data) {
  ratio <- SplitCombine.MoveRatio(params, data)
  dk1 <- ratio/(1+ratio)
  
  if (runif(1) < dk1) {
    return(AttemptMove.Split(params, data))
  } else {
    return(AttemptMove.Combine(params, data))
  }
}

GenerateBirthDeathU <- function() {
  c(runif(1), PriorDrawMu(), runif(1, min=0, max=30))
}
LogProbabilityBirthDeathU <- function(u) {
  sum(dunif(u[1], log=TRUE), LogPriorMu(u[2]), dunif(u[3], min=0, max=30, log=TRUE))
}
GenerateBirthProposal <- function(params, u) {
  ExpandParameters(params)
  
  new.k <- k+1
  new.weight <- c(weight*(1-u[1]), u[1])
  new.mu <- c(mu, circular(u[2], template="clock24"))
  new.kappa <- c(kappa, u[3])
  
  proposal <- modifyList(params, list(weight=new.weight, mu=new.mu, kappa=new.kappa, k=new.k))
  
  return(proposal)
}
#' Picks a component to remove
#' Calculates the associated helper vector $\bm u$
GenerateDeathProposal <- function(params, kj) {
  # kj is index of component to remove
  
  ExpandParameters(params)
  
  new.k <- k-1
  new.weight <- weight[-kj]
  new.weight <- new.weight / sum(new.weight) # Sum to 1
  new.mu <- mu[-kj]
  new.kappa <- kappa[-kj]
  
  new.allocation <- allocation
  new.allocation[new.allocation > kj] <- new.allocation[new.allocation > kj]-1
  
  proposal <- modifyList(params, list(weight=new.weight, mu=new.mu, kappa=new.kappa, z=new.allocation, k=new.k))
  
  return(proposal)
}
AttemptMove.BirthDeath <- function(params, data) {
  ExpandParameters(params)
  
  allocation.count <- sapply(1:k, function(j) sum(allocation == j))
  empty.components <- which(allocation.count == 0)
  
  # bk <- (length(empty.components)+1)*0.5
  bk <- 0.5
  dk1 <- (1-bk)
  move.prob.ratio <- dk1 / bk
  log.move.prob.ratio <- log(move.prob.ratio)
  
  perform.birth <- (runif(1) < dk1)

  if (perform.birth) {
    log.info("Birth ( k =", k, ")\n") # Will be appended to with MH ratio
    u <- GenerateBirthDeathU()
    proposal <- GenerateBirthProposal(params, u)
    jacobian <- (1 - u[1])^k
  } else {
    log.info("Death ( k =", k, ", empty =", length(empty.components), ")\n")
    if (length(empty.components) == 0) return(params)
    kj <- sample(empty.components, 1)
    u <- c(weight[kj], mu[kj], kappa[kj])
    proposal <- GenerateDeathProposal(params, kj)
    jacobian <- (1 - u[1])^(k-1)
  }
  
  log.u.prob <- LogProbabilityBirthDeathU(u)
  
  log.proposal <- LogPosteriorProbability(proposal, data)
  log.params <- LogPosteriorProbability(params, data)
  log.info("Proposal prob:", exp(log.proposal), "(log:", log.proposal, ")", "\nParams prob:", exp(log.params), "(log:", log.params, ")\n")
  
  if (log.proposal == -Inf) {
    print(proposal)
    cat("Impossible proposal\n")
  }
  
  if (perform.birth) {
    A <- exp(log.proposal - log.params + log.move.prob.ratio - log.u.prob + log(jacobian))
    mh.ratio <- min(1, A)
  } else {
    A <- exp(log.proposal - log.params - log.move.prob.ratio + log.u.prob - log(jacobian))
    mh.ratio <- min(1, A)
  }
  
  EvalMHRatio(mh.ratio, params, proposal)
}