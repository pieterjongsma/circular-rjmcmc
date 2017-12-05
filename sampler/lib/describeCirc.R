require(ggplot2)


# FUNCTION trigMoment  ----------------------------------------------------
# Calculates the (uncentered or centered) p-th sample trigonometric moment.
#   th:       A numeric vector containing the angles in the sample.
#             By default, this is calculated in radians.
#   p:        An integer containing the moment to be calculated. By default,
#             p = 1, so the first moment is calculated.
#   centre:   A boolean determining whether the centred moment is used. Default = FALSE.
#   rnd:      An integer giving the number of rounding digits.
# Returns:    A vector of length 4 containing the cosine and sine moments, as well as the
#             mean direction and mean resultant length.
trigMoment <- function (th, p = 1, centre = FALSE, rnd = 6) {

  # If we want the centred moments, we need the mean direction which can be obtained by
  # running this function again without centering and choosing the first moment.
  if (centre) {
    th_bar <- trigMoment(th, p = 1, centre = FALSE)[3]
  }

  # Calculate the cosine and sine p-th (centred) moments
  C_p_bar <- mean(cos(p * (th - ifelse(centre, th_bar, 0))))
  S_p_bar <- mean(sin(p * (th - ifelse(centre, th_bar, 0))))

  # Calculate mu_hat.
  mu_hat  <- meanDir(th)

  # Calculate the rho_hat
  rho_hat <- sqrt(C_p_bar^2 + S_p_bar^2 )

  return( round(
    c("Cosine Moment"         = C_p_bar,
      "Sine Moment"           = S_p_bar,
      "Mu"                    = mu_hat,
      "Rho"                   = rho_hat), rnd) )
}



# FUNCTION meanDir  ---------------------------------------------------
# Calculates the mean direction from a dataset.
#   th:       A numeric vector containing the angles in the sample, in radians.
#   na.rm:    Whether NA's will be removed.
# Returns:    A scalar, the mean direction in radians, in the range [0, 2pi].
meanDir <- function (th, na.rm=TRUE) {
  C <- sum(cos(th), na.rm=na.rm)
  S <- sum(sin(th), na.rm=na.rm)
  atan2(S, C)%%(2*pi)
}


biasCorrectionKappa <- function(kp_hat, n) {
  if (kp_hat < 2) {
    return(max(kp_hat - (2 / (n * kp_hat)), 0))
  } else {
    return(((n-1)^3 * kp_hat) / (n^3 + n))
  }
}
#
# kp_true <- .1
# n <- 30
invBiasCorrectionKappa  <- function(kp_true, n) {
  uniroot(function (x) biasCorrectionKappa(x, n) - kp_true,
        interval = c(0, 100))$root
}

# outer(c(.1, 4, 32), c(5, 30, 100), FUN = invBiasCorrectionKappa)

# Calculates a large variety of parameter estimates from a sample of angles.
#   th:       A numeric vector containing the angles in the sample.
#             By default, this is calculated in radians.
#   type:     The scale th is provided in. Can be degrees ("d") or radians ("r").
#   plot:     Whether a plotCircular() plot should be made.
# Returns:    A matrix with three moments and a matrix containing all estimates.
summaryCirc <- function (th, type = "radians", plot = FALSE) {

  if (type == "degrees" | type == "d") {
    th <- th * (pi/180)
  } else if ( ! (type == "radians" | type == "r")) {
    stop(paste("Type", type, "is unknown."))
  }

  uncenteredMomentMatrix <- rbind("p=1" = trigMoment(th, p = 1, centre = FALSE),
                                  "p=2" = trigMoment(th, p = 2, centre = FALSE),
                                  "p=3" = trigMoment(th, p = 3, centre = FALSE))
  centeredMomentMatrix   <- rbind("p=1" = trigMoment(th, p = 1, centre = TRUE),
                                  "p=2" = trigMoment(th, p = 2, centre = TRUE),
                                  "p=3" = trigMoment(th, p = 3, centre = TRUE))
  # Sample size
  n         <- length(th)

  # Mean resultant length)
  R_bar     <- uncenteredMomentMatrix[1,4]

  # Mean direction in terms of the angle it has, in radians. 2.9
  th_bar    <- uncenteredMomentMatrix[1,3]

  # Sample median direction. 2.32
  dth <- numeric(n)

  # Get the value of d(th) for each of the theta's in the sample.
  for(i in 1:n) {
    test_th <- th[i]
    dth[i]  <- pi - mean(abs(pi - abs(th - test_th)))
  }

  # The median is the value which is the minimum. Note that this assumes that
  # the data mass does not overlap 0 degrees, which it often totally does.
  th_wave   <- th[which.min(dth)]

  # Mean deviation, as in 2.33
  mean_deviation <- min(dth)


  # Sample Circular Variance. 2.12
  V         <- 1 - (R_bar)

  # Sample Circular Standard deviation. 2.13
  v         <- (-2 * log(1 - (R_bar))) ^ (1/2)

  # Second centered cosine moment
  rho_hat_2 <- centeredMomentMatrix[2 , "Cosine Moment"]

  # Delta_hat, sample circular dispersion. 2,28
  del_hat   <- (1 - rho_hat_2) / (2 * R_bar ^2)

  # Skewness, kurtosis 2.29, 2.30
  s_hat     <- (rho_hat_2 * sin(  centeredMomentMatrix[2, "Mu"] - 2 * centeredMomentMatrix[1, "Mu"])) /
    ((1 - R_bar) ^ (3/2))
  k_hat     <- (rho_hat_2 * cos(uncenteredMomentMatrix[2, "Mu"] - 2 * centeredMomentMatrix[1, "Mu"])- (R_bar^4) ) /
    ((1 - R_bar) ^ (2))


  estMatrix <- rbind(
    "Mean direction"      = th_bar,
    "Median direction"    = th_wave,
    "Mean deviation"      = mean_deviation,
    "Sample Circular Var" = V,
    "Sample Circular SD"  = v,
    "Dispersion"          = del_hat,
    "Skewness"            = s_hat,
    "Kurtosis"            = k_hat

  )

  colnames(estMatrix) <- "Value"

  if (plot) plotCircular(th, plotMoments = 3)

  return(
    list (
      "Uncentered sample moments" = uncenteredMomentMatrix,
      "Centered sample moments"   = centeredMomentMatrix,
      Estimates = estMatrix
    )
  )
}









###
### FUNCTION plotCircular
###
#
#   th:           A numeric vector containing the angles in the sample.
#                 By default, this is calculated in radians.
#   alpha:        A numeric determining opaqueness of plotted points.
#   type:         A string that is either "radians" or "degrees", the scale of th.
#   plotMean:     A boolean determining whether the mean direction should be plotted.
#   plotMoments:  A numeric determining how many moments should be plotted.
# Returns:    A vector of length 4 containing the cosine and sine moments, as well as the
#             mean direction and mean resultant length.
plotCircular <- function(th, alpha = 0.35, inputtype = "radians",
                         outputtype = "radians", plotMean = TRUE, plotMoments = 0, col = rgb(0,0,0, alpha)) {

  if (!require(plotrix)) stop("\n Package 'plotrix' must be installed! \n")

  if (inputtype == "degrees") th <- th*(pi/180)

  th_bar <- meanDir(th)

  # Empty plot
  plot(-1:1, -1:1, type = "n",
       xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
       xlab = "", ylab = "",
       asp = 1, xaxt = "n", yaxt = 'n', bty = 'n')

  draw.circle(0,0, radius = 1.05)

  segments(0,0, 1.05, 0, col = "gray")

  # Plot the points
  points(cbind(cos(th)), sin(th), pch = 16, cex = 1.3, col=col)

  # Plot the mean
  if (plotMean) {

    points(cos(th_bar), sin(th_bar), col = "skyblue3", pch = 19)
    segments(0,0, cos(th_bar), sin(th_bar), col = "skyblue3")

    text(0.6, 0.07, labels = paste("0", outputtype), cex = 0.6, col = "gray")
    text(0, -1.15,
         labels = paste("Mean direction =",
                        round(ifelse(outputtype == 'radians', th_bar, th_bar*(180/pi)), 3),
                        outputtype),
         col = "skyblue3", cex = 0.8)
  }

  # Plot moments
  if (plotMoments > 0) {
    for (moment in 1:plotMoments) {
      coord <- trigMoment(th, p = moment, centre = FALSE)[1:2]
      text(coord[1], coord[2], moment, cex = .6, col = "dark green")
    }
  }
}




###
### FUNCTION ggplotCircular
###
#
#   th:           A numeric vector containing the angles in the sample.
#                 By default, this is calculated in radians.
#   alpha:        A numeric determining opaqueness of plotted points.
#   type:         A string that is either "radians" or "degrees", the scale of th.
#   plotMean:     A boolean determining whether the mean direction should be plotted.
#   plotMoments:  A numeric determining how many moments should be plotted.
# Returns:    A vector of length 4 containing the cosine and sine moments, as well as the
#             mean direction and mean resultant length.
ggplotCircular <- function(th, alpha = 0.35, inputtype = "radians",
                         outputtype = "radians", plotMean = TRUE,
                         col = rgb(0,0,0, alpha)) {

  if (!require(plotrix)) stop("\n Package 'plotrix' must be installed! \n")

  if (inputtype == "degrees") th <- th*(pi/180)

  th_bar <- meanDir(th)

  csq <- seq(0, 2*pi, length.out = 180)
  xsq <- cos(csq)
  ysq <- sin(csq)
  dsq <- data.frame(xsq, ysq)

  p <- ggplot(aes(x=xsq, y=ysq), data=dsq) +
    geom_path() + coord_fixed() +
    theme(panel.background = element_blank(),
          line = element_blank(),
          text = element_blank()) +
    geom_path(aes(x=x, y=y),
              data=data.frame(x=c(0, 1.05), y=c(0, 0)), colour="gray") +
    geom_point(aes(x=cos(th), y=sin(th)), size=4, alpha=alpha, data=data.frame(th))

  # Plot the mean
  if (plotMean) {
    meanTxt <- paste("Mean direction:",
                     round(ifelse(outputtype == 'radians',
                                  th_bar,
                                  th_bar*(180/pi)), 3), outputtype)

    p  <- p + geom_point(aes(cos(th_bar), sin(th_bar)),
                         data=data.frame(th_bar), color= "skyblue3", pch = 19) +
      geom_path(aes(x=x, y=y),
                data=data.frame(x=c(0, cos(th_bar)), y=c(0, sin(th_bar))),
                color= "skyblue3") +
      geom_text(aes(x=x, y=y), label=meanTxt,
                data=data.frame(x=0, y=-1.15), color= "skyblue3", size=5)
  }
  p
}





# Adds a vector of angels to the plot.
addAngle <- function(th, col="tomato") {
  for (i in 1:length(th)) {
    points(cos(th[i]), sin(th[i]), col = col, pch = 19)
    segments(0,0, cos(th[i]), sin(th[i]), col = col)
  }
}



###
### FUNCTION plotCircularDensity
### Plots a circular density on a circle.
###
plotCircularDensity <- function(FUN, r=.6, res=360, add=FALSE, col="black", ...) {

  if (!require(plotrix)) stop("\n Package 'plotrix' must be installed! \n")

  require(plotrix)
  if (!add) {
    plot(-1:1, -1:1, type = "n",
         xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2),
         xlab = "", ylab = "",
         asp = 1, xaxt = "n", yaxt = 'n', bty = 'n')
  }

  draw.circle(0,0, radius = r)
  segments(0,0, r, 0, col = "gray")

  pts <- seq(from=0, to=360, by=360/res)*(pi/180)
  d <- 1+sapply(pts, FUN)

  cdp <- r*d*cos(pts)
  sdp <- r*d*sin(pts)

  lines(cdp, sdp, col=col)
}






###
### FUNCTION circularQuantile
###
# A wrapper for quantile() that first reverse-centres the circular data in order
# to prevent the results being influenced by an arbitrary starting point.
#   th:         A numeric vector containing the angles in the sample.
#   ...:        Arguments to be passed to quantile().
# Returns:    A vector containing the angles of the desired quantiles.
circularQuantile <- function (th, ...) {

  # The desired rotation before taking the quantile.
  rotation <- pi-meanDir(th)

  # Centre the data, and move it as far away from 0 radians as possible by
  # th+rot. Then, apply the quantile function, and rotate back.
  quantile(th+rotation, ...)-rotation
}




# Very similar to circularquantile, tries to replace a set of angles with values
# above 0 and below 2pi to a set with values above 0 and below 0, for more
# pretty plotting, for instance.
setCircularRangeContinuous <- function (th) {

  # The desired rotation.
  rotation <- pi-meanDir(th)

  # Centre the data, and move it as far away from 0 radians as possible by
  # th+rot. Then, apply the quantile function, and rotate back.
  (th+rotation)%%(2*pi) - rotation
}







