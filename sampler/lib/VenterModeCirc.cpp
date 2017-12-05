/*
----------------------------------------------------------
VenterMode.cpp
Calculate the Highest Posterior density (HPD, as in Venter (1967), to either
estimate the mode or obtain the shortest credible interval, for example.

Kees Tim Mulder
Last updated: November 2014

This work was supported by a Vidi grant awarded to I. Klugkist from the
Dutch Organization for Scientific research (NWO 452-12-010).
----------------------------------------------------------
*/


#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
double hmodecirc(NumericVector x, double cip) {
  /* FUNCTION hmode -------------------------------------------
  Estimate the mode by finding the highest posterior density interval.

  x:      Sample from which to estimate the mode.
  cip:    Bandwith for the algorithm, ranging from 0 to 1.

  Returns: A scalar containing the estimate of the mode.
  ------------------------------------------------------------ */

  int n, cil, chiv;
  double ln, M;

  n = x.size();
  NumericVector sx = clone(x);
  std::sort(sx.begin(), sx.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Size of the currently smallest interval.
  ln = sx[cil]-sx[0];

  for (int i=0; i < (n-cil); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    double ln_new = (sx[(i+cil) % n]-sx[i]);
    if (ln_new < 0) ln_new = 2*M_PI - sx[i] + sx[(i+cil) % n];
    if (ln > ln_new) {
      ln = ln_new;
      chiv = i;
    }
  }

  M = (sx[chiv+cil]+sx[chiv])/2;

  return M;
}


// [[Rcpp::export]]
NumericVector hmodecicirc(NumericVector x, double cip) {
  /* FUNCTION hmodecicirc -------------------------------------
  Find the highest posterior density interval.

  x:      Sample from which to estimate the interval.
  cip:    Bandwith for the algorithm, ranging from 0 to 1.

  Returns: An vector of length 2 containing
           lower and upper bound of the interval.
  ------------------------------------------------------------ */

  int n, cil, chiv;
  double ln;

  n = x.size();
  NumericVector sx = clone(x);
  std::sort(sx.begin(), sx.end());

  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Length of the currently smallest interval.
  ln = sx[cil]-sx[0];

  for (int i=0; i < n; i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    double ln_new = (sx[(i+cil) % n]-sx[i]);
    if (ln_new < 0) ln_new = 2*M_PI - sx[i] + sx[(i+cil) % n];
    if (ln > ln_new) {
      ln = ln_new;
      chiv = i;
    }
  }

  NumericVector M(2);
  M[0] = sx[chiv];
  M[1] = sx[(chiv+cil) % n];

  return M;
}
