#include <Rcpp.h>
using namespace Rcpp;

IntegerVector rmultinom_1(unsigned int &size, NumericVector &probs, unsigned int &N) {
  IntegerVector outcome(N);
  rmultinom(size, probs.begin(), N, outcome.begin());
  return outcome;
}

// [[Rcpp::export]]
IntegerMatrix rmultinom_rcpp(unsigned int &n, unsigned int &size, NumericVector &probs) {
  unsigned int N = probs.length();
  IntegerMatrix sim(N, n);
  for (unsigned int i = 0; i < n; i++) {
    sim(_,i) = rmultinom_1(size, probs, N);
  }
  return sim;
}