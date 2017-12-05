/*
----------------------------------------------------------
VMMHC.cpp
This is an RCPP implementation the  of the algorithm presented by Forbes & Mardia (2014), adapted to
sampling from the posterior of multiple groups of data from the von Mises distribution.

Kees Tim Mulder
Last updated: November 2014

This work was supported by a Vidi grant awarded to I. Klugkist from the
Dutch Organization for Scientific research (NWO 452-12-010).
----------------------------------------------------------
*/

#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>

// [[Rcpp::export]]
NumericVector rvmc(int n, double mu, double kp) {
  /* FUNCTION rvmc -------------------------------------------
  Generate random variates from the von Mises distribution.

  n:      The number of random variates required.
  mu:     The required mean direction, mu.
  kp:     The required concentration, kappa.

  Returns: A vector of length n containing VM random variates.
  ------------------------------------------------------------ */

  // If kappa is very small, return a circular uniform draw, as otherwise the
  // algorithm will fail.
  if (kp < .0000001) {
      return runif(n, 0, 4.0*atan(1));
  }

  // Main algorithm.
  NumericVector th(n);
  int sn;
  double a, b, r, u1, u2, u3, z, f, c;
  bool cont;

  // Setup.
  a = 1 + sqrt(1 + 4.0 * pow(kp, 2));
  b = (a - (sqrt(2.0*a)))/(2.0*kp);
  r = (1 + pow(b,2))/(2.0*b);

  for (int i=0; i<n; i++) {

    cont = TRUE;

    do {
      u1 = runif(1, 0, 1)[0];
      u2 = runif(1, 0, 1)[0];
      u3 = runif(1, 0, 1)[0];

      // STEP 1
      z = cos(4*atan(1)*u1);
      f = (1 + r*z)/(r + z);
      c = kp * (r - f);

      // STEP 2
      if (c*(2-c) - u2 > 0) cont=FALSE;

      // STEP 3
      if (log(c/u2) + 1 - c >= 0) cont=FALSE;
    } while (cont);

    // STEP 4
    if (u3 - 0.5 > 0) {
      sn = 1;
    } else {
      sn = -1;
    }

    th[i] = fmod(sn * acos(f) + mu, 8.0*atan(1));
}

  return th;
}


// Approximation of Lamberts W.
double Wapprox (double t) {
  return exp(1)*t / (1 + 1 / (pow(2 * exp(1) * t + 2, -0.5) + 1 / (exp(1) - 1) - pow(2, -0.5)));
}

// The rejection sampler
// [[Rcpp::export]]
Rcpp::List FMC(Rcpp::List th, double kp_start,
               NumericVector mu_n, NumericVector R_n, double R_t, int m_t,
               int Qb, int lag) {
  /* FUNCTION FMC -------------------------------------------
  Generates samples from the posterior of k von Mises distributions, each with
  a separate mean, but with one common concentration kappa.
  ------------------------------------------------------------ */

  int k = mu_n.size();

  NumericMatrix mu(Qb, k);
  NumericVector kp(Qb);

  // The variable that saves the current value of mu, which
  // is only copied if the current iteration is not thinned
  // out.
  NumericVector mu_new(k);

  int eta = m_t;
  int Qbylag = Qb * lag;
  int idvlag = 1;
  int candidates = 0;

  double etab0, b0, kl, ku, c1, c2, c3, c4, c5, c6, k0, i0, r, beta, eps,
  alph, x, kp_new, kp_cur, u, v1, v2, v;

  // Indicating whether the candidate is accepted.
  bool cont;

  kp_new = kp_start;

  for (int i = 0; i < Qbylag; i++) {
    kp_cur = kp_new;

    // Analytic mean for mu
    for (int j = 0; j < k; j++) {
      mu_new[j] = rvmc(1, mu_n[j], R_n[j]*kp_cur)[0];
    }

    // Setup: Calculate beta_0, here called b0. This is also called beta_t in my
    // paper, which represents the fact that it is bundled from different groups
    // and includes prior information.
    etab0 = 0;
    for (int j = 0; j < k; j++) {
      etab0 -= R_n[j] * cos(mu_new[j] - mu_n[j]);
    }
    b0 = etab0 / eta;

    // Setup: Compute approximately optimal parameters for the rejection
    // sampler. kappa_0 is called k0.
    kl   = 2 / (etab0 + sqrt(2 * eta + etab0 * etab0));
    ku   = (2 + (1 / eta) ) / (etab0 + b0 + sqrt(2 * eta + 1 + etab0 * etab0));
    c1   = 0.5 + 0.5 * (1 - 0.5 / eta) / eta;
    k0   = (1 - c1) * kl + c1 * ku;
    i0   = boost::math::cyl_bessel_i(0, k0);
    r    = boost::math::cyl_bessel_i(1, k0) / i0;
    c2   = 0.25 / eta - 2 / (3 * sqrt(eta));
    if (b0 < c2) {
      beta = b0 + 1;
    } else {
      beta = b0 + r + (1 - r) / (1 + 40 * eta * pow(b0 - c2, 2));
    }
    c3   = (log(i0) / k0 - beta + b0) / (beta - b0 - r);
    c4   = Wapprox(c3 * exp(c3));
    eps  = c4 * k0 / (c3 - c4);
    alph = (beta - b0 - r) * (k0 + eps);
    c5   = log(i0);

    // Apply rejection sampler
    cont = TRUE;
    do {

      // Save the number of candidates
      candidates = candidates + 1;

      // Draw values from the gamma distribution with the
      // tweaked parameters
      x = rgamma(1, eta * alph + 1, 1.0 / (eta * beta))[0];

      if (x > eps) {

        // Compute kp_new and v
        kp_new = x - eps;
        c6 = 0.5 * log(8 * atan(1) * kp_new) - kp_new;
        u  = runif(1, 0, 1)[0];
        v1 = log(u) / eta - (beta - b0) * (kp_new - k0);
        v2 = alph * log((kp_new+eps) / (k0 + eps)) - c5;
        v  = v1 + v2;

        // Break the loop if these tests are passed.
        if (kp_new < 0.258 || v < c6) {
          if ( v < c6 - log(1 + 1 / 2 * kp_new) ||
              v < -log(boost::math::cyl_bessel_i(0, kp_new))) {
            cont = FALSE;
          }
        }
      }
    } while (cont);

    //    For non-thinned out iterations, save the current value.
    if (i % lag == 0) {
      idvlag = i/lag;
      mu(idvlag, _) = mu_new;
      kp[idvlag] = kp_new;
    }
  }

  //  Gather output.
  NumericMatrix out(Qb, 1 + k);
  out( _, 0) = kp;
  for (int mui = 0; mui < k; mui++) {
    out( _, mui + 1) = mu( _, mui);
  }

  return Rcpp::List::create(Rcpp::Named("sam") = out,
         Rcpp::Named("att") = Rcpp::wrap(candidates));
}

