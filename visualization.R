require('circular')

source('configuration.R')
source('sampler/circular_sampler.R')

PlotWeight <- function(weight, color.offset=0) {
  iterations <- nrow(weight)
  
  plot(0, 0, type="n", xlim=c(1, iterations), ylim=c(0, 1),
    ylab=expression(w), xlab="Iteration")
  for (j in 1:ncol(weight)) {
    weights <- weight[, j]
    lines(1:iterations, weights, col=custom.colors[[j+color.offset]])
  }
}

PlotMus <- function(mu) {
  par(mfrow=c(ncol(mu), 1), mar=c(3, 3, 3, 3))
  
  for (j in 1:ncol(mu)) {
    PlotMu(as.matrix(mu[, j]), color.offset=j-1)
  }
}
PlotMu <- function(mu, color.offset=0) {
  iterations <- nrow(mu)
  
  plot(0, 0, type="n", xlim=c(1, iterations), ylim=c(-pi, pi),
    ylab=expression(mu), xlab="Iteration")
  for (j in 1:ncol(mu)) {
    mus <- mu[, j]
    lines(1:iterations, mus, col=custom.colors[[j+color.offset]])
  }
}
NormalizeMu <- function(mu) {
  mu[mu < mean(mu)-pi] <- mu[mu < mean(mu)-pi] + 2*pi
  mu[mu > mean(mu)+pi] <- mu[mu > mean(mu)+pi] - 2*pi
  mu
}

PlotKappa <- function(kappa, color.offset=0) {
  iterations <- nrow(kappa)
  
  plot(0, 0, type="n", xlim=c(1, iterations), ylim=c(0, 30),
    ylab=expression(kappa), xlab="Iteration")
  for (j in 1:ncol(kappa)) {
    kappas <- kappa[, j]
    lines(1:iterations, kappas, col=custom.colors[[j+color.offset]])
  }
}


Hist <- function(data, ...) {
  data.joined <- as.vector(data)
  hist(data.joined, breaks=20, ...)
}
HistWeight <- function(weight) { Hist(weight, xlim=c(0, 1)) }
HistMu <- function(mu) { Hist(mu, xlim=c(0, 2*pi)) }
HistKappa <- function(kappa) { Hist(kappa, xlim=c(0, 30)) }


CurveForParameters <- function(params, ...) {
  ExpandParameters(params)
  
  dens <- function(xs) {
    sapply(xs, function(x) {
      sum(sapply(1:k, function(j) {
        weight[j] * dvonmises(x, mu=mu[j], kappa=kappa[j])
      }))
    })
  }
  
  curve.circular(dens, ...)
}

CurvesForParameters <- function(params, range=NA, ...) {
  ExpandParameters(params)
  
  if (is.na(range)) range <- 1:k
  
  for (j in range) {
    dens <- function(xs) {
      sapply(xs, function(x) {
        weight[j] * dvonmises(x, mu=mu[j], kappa=kappa[j])
      })
    }
  
    curve.circular(dens, ...)
  }
}



custom.colors <- list(
  turquoise="#1abc9c",
  green.sea="#16a085",
  sun.flower="#f1c40f",
  orange="#f39c12",
  emerald="#2ecc71",
  nephritis="#27ae60",
  carrot="#e67e22",
  pumpkin="#d35400",
  peter.river="#3498db",
  belize.hole="#2980b9",
  alizarin="#e74c3c",
  pomegranate="#c0392b",
  amethyst="#9b59b6",
  wisteria="#8e44ad",
  clouds="#ecf0f1",
  silver="#bdc3c7",
  wet.asphalt="#34495e",
  midnight.blue="#2c3e50",
  concrete="#95a5a6",
  asbestos="#7f8c8d"
)