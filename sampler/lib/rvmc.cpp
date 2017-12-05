#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rvmc(int n, double mu, double kp) {
  //
  //    FUNCTION rvmc
  //
  // This is an RCPP implementation of the algorithm for generating Von Mises
  // random variates on page 49 of N.I. Fishers' Statistical Analysis of
  // Circular data.
  //
  //  n:      The number of random variates required.
  //  mu:     The required mean direction, mu.
  //  kp:     The required concentration, kappa.
  //
  // Returns: A vector of length n containing VM random variates.

  // If kappa is very small, return a circular uniform draw, as otherwise the
  // algorithm will fail.
  if (kp < .0000001) {
      return runif(n, 0, 8.0*atan(1));
  }

  NumericVector th(n);
  int sn;
  double a, b, r, u1, u2, u3, z, f, c;
  bool cont;

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

      // STEP 3 Smaller than is now larger eq than, is turned around, as in the
      // loop going back to step 1 is automatic, while if not, we must break and
      // go to step 4.
      if (log(c/u2) + 1 - c >= 0) cont=FALSE;
    } while (cont);

    //# STEP 4
    if (u3 - 0.5 > 0) {
      sn = 1;
    } else {
      sn = -1;
    }

    th[i] = fmod(sn * acos(f) + mu, 8.0*atan(1));
}

  return th;
}
