# FUNCTION dvm  -----------------------------------------------------------
# The probability density function of the von Mises distribution.
#   th:       Theta, the angle in radians for which the density will be found.
#   mu, kappa:  Location and concentration parameters of von Mises, mu and kappa.
# Returns:    The density at th.
dvm <- function (th, mu = 0, kappa = 1) {

  # Fisher, 3.3.1
  (1 / (2*pi*besselI(kappa, 0)) ) * exp (kappa * cos(th - mu))

}



# FUNCTION pvm  -----------------------------------------------------------
# The cumulative distribution function of the von Mises distribution.
# This simply uses the integrate function in R, as the cdf of von Mises is not
# analytically defined, see also eq. 3.50.
#   th:       Theta, the angle in radians for which the density will be found.
#   kappa:      Concentration parameter of von Mises, kappa.
# Returns:    The cumulative probability at th.
pvm <- function (th, kappa = 1) {
  integrate(dvm, -pi, th, kappa = kappa)$value
}



# FUNCTION inv_F_k  -------------------------------------------------------
# The inverse cumulative distribution function of the von Mises distribution.
# This is calculated by an algorithm in Fisher.
#   P:        The cumulative probability for which the quantile will be found.
#   kappa:      Concentration parameter of von Mises, kappa.
# Returns:    The quantile in radians that corresponds to the cumulative
#             probability P.
inv_F_k <- function (P, kappa = 1, tolerance = .00001) {

  # Page 53, in Fisher, section 3.3.6.
  ### STEP 1 ###
  f <- 0.5
  t <- 0
  c <- log(besselI(kappa, 0))

  ### STEP 2.1 ###
  repeat {
    g <- f - P
    d <- exp(log(abs(g)) + c - kappa * cos(t))

    ### STEP 3 ###
    if (d < tolerance) break

    ### STEP 2.2 ###
    t <- t - sign(g) * d
    f <- pvm(t, kappa = kappa)
  }

  ### STEP 4 ###
  t
}



# FUNCTION qvm  -------------------------------------------------------
# Calculates the theoretical quantiles of the von Mises distribution for a
# vector of probabilities.
#   th:       Theta, the angle in radians for which the density will be found.
#   kappa:      Concentration parameter of von Mises, kappa.
# Returns:    The quantiles in radians that corresponds to the cumulative
#             probability P.
qvm <- function (p, kappa = 1) {
  q <- numeric(length(p))
  for (i in 1:length(p)) q[i] <- inv_F_k(p[i], kappa = kappa)
  q
}



# FUNCTION qvm  -------------------------------------------------------
# Calculate a set of von Mises quantiles from a sample size.
#   n:        The sample size for which to obtain quantiles.
#   kappa:      Concentration parameter of von Mises, kappa.
# Returns:    The quantiles in radians that corresponds to the cumulative
#             probability P.
qvm_n <- function (n, kappa = 1) {
  qvm(1:n/(n+1), kappa = kappa)
}



# FUNCTION approxKappaML  -----------------------------------------------
# Approximate kappa in a sample from the von Mises distribution.
#   th:       The vectors in the sample.
# Returns:    An approximation for the maximum likelihood estimate of kappa.
approxKappaML <- function (th) {
  R_bar <- sqrt(mean(cos(th))^2 + mean(sin(th))^2)

  n <- length(th)

  kappaML <- 0

  # 4.40, 4.41
  if (R_bar < .53) {
    kappaML <- 2 * R_bar + R_bar ^ 3 + 5/6 * R_bar ^ 5
  } else if (R_bar < .85) {
    kappaML <- -.4 + 1.39 * R_bar + .43/(1-R_bar)
  } else {
    kappaML <- 1/(R_bar^3 - 4 * R_bar^2 + 3 * R_bar)
  }

  if (n <= 15 & kappaML < 2){
    kappaML <- max(kappaML - (1/(2*n*kappaML)), 0)
  } else if (n <= 15 & kappaML >= 2) {
    kappaML <- ((n-1)^3 * kappaML)/(n^3 + n)
  }
  kappaML
}








# # FUNCTION kappaML  -----------------------------------------------
# # Approximate kappa in a sample from the von Mises distribution.
# #   th:       The vectors in the sample.
# # Returns:    An approximation for kappa.
# estKappaML <- function (th, tol = .01) {
#   besselA1 <- function(x) {
#     besselI(x, 1) / besselI(x, 0)
#   }
#
#
#   R_bar <- sqrt(mean(cos(th))^2 + mean(sin(th))^2)
#
#   #   besselA1(kappaML) == R_bar
#
#   kappaML <- approxKappaML(th)
#
#   repeat {
#
#     if ((besselA1(kappaML) - R_bar) < tol) break
#
#
#   }
#
# }





# CHECKS ------------------------------------------------------------------

#
# setwd("C:/Dropbox/Master's Thesis/")
# source("Code/General/sampleCirc.R")
# source("Code/General/describeCirc.R")
# require(circular)
#
# # Checks the pdf
# plot(dvm, xlim = c(-pi, pi), xlab = "th", ylab = "f(th)")
#
# # Checks the cdf
# th <- (0:620/100)-3.1
# plot(th, sapply(th, pvm),
#      xlim = c(-pi, pi),
#      ylim = c(0, 1.1),
#      type = "l",
#      xlab = "Angles (Radians)",
#      ylab = "Cumulative distribution function")
#
#
# p <- c(0, 1:49/50, 1)
# # Checks the inverse cdf
# plot(p, sapply(p, inv_F_k, kappa = 1),
#      xlim = c(0, 1),
#      ylim = c(-pi, pi),
#      type = "l",
#      xlab = "Cumulative probability",
#      ylab = "Angles (Radians)")
#
# th <- rvm(100, mu=0, kappa = 5)
# th2 <- as.numeric(rvonmises(100, mu = 0, kappa = 5))
#
# z <-  sort(sin(0.5 * (th)))
# z2 <-  sort(sin(th2))
#
#
# vmquantiles <- sin(0.5 * qvm_n(length(th), kappa = 5))
#
#
# n <- 100
# vmq <- as.numeric(qvonmises(1:n/(n+1), mu = 0, kappa = 5))
# vmquantiles2 <- sin(vmq)
#
#
# # plot(rvonmises(100, mu = 0, kappa = 5))
# # plot(circular(th))
#
#
# plot(vmquantiles2, z2,
#      xlim = c(-1, 1),
#      ylim = c(-1, 1),
#      type = "p",
#      xlab = "von Mises quantiles",
#      ylab = "Sample Quantiles")
# abline(0,1)
#
#
# summaryCirc(th, plot = TRUE)
#
# # Check if we can use the ML of kappa
# th <- rvm(10000, mu = 2, kappa = 2.6)
# summaryCirc(th)
# besselI(2.6, 1)/besselI(2.6, 0)
# # Dit is gelijk aan de R_bar (rho) van p=1, en dat klopt.
