#include <Rcpp.h>
using namespace Rcpp;

IntegerVector rmultinom_1(unsigned int &size, NumericVector &probs, unsigned int &N) {
  IntegerVector outcome(N);
  rmultinom(size, probs.begin(), N, outcome.begin());
  return outcome;
}

// [[Rcpp::export]]
IntegerMatrix rmultinomSizeVecProbMat_rcpp(IntegerVector &size, NumericMatrix &probs) {
  unsigned int N = probs.nrow();
  unsigned int n = probs.ncol();
  IntegerMatrix sim(N, n);
  NumericVector theseProbs(N);
  unsigned int thisSize;
  for (unsigned int i = 0; i < n; i++) {
    theseProbs = probs(_,i);
    thisSize = size[i];
    sim(_,i) = rmultinom_1(thisSize, theseProbs, N);
  }
  return sim;
}

// [[Rcpp::export]]
IntegerMatrix rmultinomSizeVec_rcpp(IntegerVector &size, NumericVector &probs) {
  unsigned int N = probs.length();
  unsigned int n = size.length();
  IntegerMatrix sim(N, n);
  unsigned int thisSize;
  for (unsigned int i = 0; i < n; i++) {
    thisSize = size[i];
    sim(_,i) = rmultinom_1(thisSize, probs, N);
  }
  return sim;
}