// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>

#include <Rcpp.h>
#include <math.h>

using namespace Rcpp;

const double pi = boost::math::constants::pi<double>();




// [[Rcpp::export]]
double Wapprox (double t) {
  // Approximation of Lamberts W.
  return exp(1)*t /
  (1 + 1 / (pow(2 * exp(1) * t + 2, -0.5) + 1 / (exp(1) - 1) - pow(2, -0.5)));
}

// [[Rcpp::export]]
NumericVector sampleKappa(double etag, int eta) {
  // Function to sample values for kappa in a von Mises distribution. Etag should
  // be -R * cos(mu - theta_bar). eta is the posterior n, which is n + c where c is
  // the number of observations contained in the conjugate prior. For
  // uninformative, c = 0 and eta = n.

  // beta_0 in Forbes & Mardia (2014) is renamed g here to avoid confusion with
  // the intercept in the GLM model.
  double g, kl, ku, c1, c2, c3, c4, c5, c6, beta, eps,
  alph, x, u, v1, v2, v;

  long double k0, i0, r, kp_can;
  
  int nAttempts = 0;

  try {
    // Boolean indicating whether the current candidate is accepted.
    bool cont;

    g = etag / eta;

    // Setup: Compute approximately optimal parameters for the rejection
    // sampler. kappa_0 is called k0.
    kl   = 2 / (etag + sqrt(2 * eta + etag * etag));
    ku   = (2 + (1 / eta) ) / (etag + g + sqrt(2 * eta + 1 + etag * etag));
    c1   = 0.5 + 0.5 * (1 - 0.5 / eta) / eta;
    k0   = (1 - c1) * kl + c1 * ku;
    i0   = boost::math::cyl_bessel_i(0, k0);
    r    = boost::math::cyl_bessel_i(1, k0) / i0;
    c2   = 0.25 / eta - 2 / (3 * sqrt(eta));
    if (g < c2) {
      beta = g + 1;
    } else {
      beta = g + r + (1 - r) / (1 + 40 * eta * pow(g - c2, 2));
    }
    c3   = (log(i0) / k0 - beta + g) / (beta - g - r);
    c4   = Wapprox(c3 * exp(c3));
    eps  = c4 * k0 / (c3 - c4);
    alph = (beta - g - r) * (k0 + eps);
    c5   = log(i0);

    // Apply rejection sampler.
    cont = TRUE;
    do {

      // Draw values from the gamma distribution with the
      // tweaked parameters.
      x = rgamma(1, eta * alph + 1, 1.0 / (eta * beta))[0];

      // Save the number of candidates.
      nAttempts++;

      if (x > eps) {

        // Compute kp_can and v.
        kp_can = x - eps;
        c6 = 0.5 * log(2 * pi * kp_can) - kp_can;
        u  = runif(1, 0, 1)[0];
        v1 = log(u) / eta - (beta - g) * (kp_can - k0);
        v2 = alph * log((kp_can+eps) / (k0 + eps)) - c5;
        v  = v1 + v2;

        // Break the loop if these tests are passed.
        if (kp_can < 0.258 || v < c6) {
          if ( v < c6 - log(1 + 1 / 2 * kp_can) ||
               v < -log(boost::math::cyl_bessel_i(0, kp_can))) {
            cont = FALSE;
          }
        }

      }
    
      if (nAttempts > 20) {
        kp_can = -1;
        cont = FALSE;
      }
    } while (cont);
  } catch (...) {
    kp_can = -1;
  }

  NumericVector out = NumericVector(2);
  out(0) = kp_can;
  out(1) = nAttempts;
  return out;
}
