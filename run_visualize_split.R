#! /usr/local/bin/Rscript --vanilla

source('configuration.R')
source('visualization.R')

frame <- 1.3

definition <- list(
  k=1,
  mu=circular(1, template="clock24"),
  kappa=0.5,
  weight=1
)
weight <- definition$weight
rho <- A1(definition$kappa)
alpha <- rho*cos(definition$mu)
beta <- rho*sin(definition$mu)

file <- "figures/original.pdf"
pdf(file, width=4, height=5)
par(mar=c(0, 0, 0, 0))
plot(circular(NA, rotation="clock", zero=pi/2), xlim=c(-frame, frame), ylim=c(-frame, frame))
CurveForParameters(definition, add=TRUE)
lines(c(0, beta), c(0, alpha))
dev.off()


weight1 <- 0.3
weight2 <- weight - weight1

alpha1 <- 0.6
beta1 <- -0.1
rho1 <- sqrt(alpha1^2 + beta1^2)

alpha2 <- (alpha*weight - alpha1*weight1)/weight2
beta2 <- (beta*weight - beta1*weight1)/weight2
rho2 <- sqrt(alpha2^2 + beta2^2)

definition1 <- list(
  k=2,
  mu=circular(c(atan2(beta1, alpha1), atan2(beta2, alpha2)), template="clock24"),
  kappa=c(A1inv(rho1), A1inv(rho2)),
  weight=c(weight1, weight2)
)

file <- "figures/split.pdf"
pdf(file, width=4, height=5)
par(mar=c(0, 0, 0, 0))
plot(circular(NA, rotation="clock", zero=pi/2), xlim=c(-frame, frame), ylim=c(-frame, frame))
CurvesForParameters(definition1, add=TRUE)
lines(c(0, beta1), c(0, alpha1))
lines(c(0, beta2), c(0, alpha2))
dev.off()
