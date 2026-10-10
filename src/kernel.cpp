#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// acc holds the pairs i < j in the order of a dist object: column i, rows j > i.
// Rows must be grouped by taxon; sample ids are 0-based.
// [[Rcpp::export]]
void accumulate_pairs(NumericVector acc, IntegerVector sample, IntegerVector taxon,
                      NumericVector quant, int n, bool product = false) {
  const R_xlen_t N = sample.size();
  const R_xlen_t size = (R_xlen_t) n * (n - 1) / 2;
  if (acc.size() != size) stop("`acc` must have n * (n - 1) / 2 elements");
  double* out = acc.begin();
  R_xlen_t start = 0;
  while (start < N) {
    R_xlen_t end = start + 1;
    while (end < N && taxon[end] == taxon[start]) ++end;
    for (R_xlen_t a = start; a < end; ++a) {
      const int sa = sample[a];
      if (sa < 0 || sa >= n) stop("`sample` ids must be in [0, n)");
      const double qa = quant[a];
      for (R_xlen_t b = a + 1; b < end; ++b) {
        const int sb = sample[b];
        if (sb < 0 || sb >= n) stop("`sample` ids must be in [0, n)");
        const int i = sa < sb ? sa : sb;
        const int j = sa < sb ? sb : sa;
        const double qb = quant[b];
        out[(R_xlen_t) i * (2 * (R_xlen_t) n - i - 1) / 2 + (j - i - 1)] +=
          product ? qa * qb : (qa < qb ? qa : qb);
      }
    }
    start = end;
  }
}

// method: 0 = bray, 1 = jaccard, 2 = euclidean. acc holds the shared minima
// (or products for euclidean) and is replaced by the distances.
// [[Rcpp::export]]
void finalize_dist(NumericVector acc, NumericVector total, NumericVector square,
                   int n, int method) {
  if (acc.size() != (R_xlen_t) n * (n - 1) / 2) stop("`acc` must have n * (n - 1) / 2 elements");
  double* out = acc.begin();
  R_xlen_t k = 0;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j, ++k) {
      if (method == 2) {
        const double d2 = square[i] + square[j] - 2 * out[k];
        out[k] = d2 > 0 ? std::sqrt(d2) : 0.0;
      } else {
        const double b = 1 - 2 * out[k] / (total[i] + total[j]);
        out[k] = method == 1 ? 2 * b / (1 + b) : b;
      }
    }
  }
}
