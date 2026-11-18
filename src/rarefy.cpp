#include <Rcpp.h>
using namespace Rcpp;

// Rows must be grouped by sample; start holds the first row of each group and
// ends with the number of rows. A group is drawn taxon by taxon, each from the
// hypergeometric distribution of what is left, which is a draw without
// replacement of `depth` individuals.
// [[Rcpp::export]]
NumericVector rarefy_runs(NumericVector quant, IntegerVector start, double depth) {
  NumericVector out = clone(quant);
  for (R_xlen_t g = 0; g + 1 < start.size(); g++) {
    double total = 0;
    for (int i = start[g]; i < start[g + 1]; i++) total += quant[i];
    if (total <= depth) continue;
    double left = depth;
    for (int i = start[g]; i < start[g + 1]; i++) {
      total -= quant[i];
      out[i] = left > 0 ? R::rhyper(quant[i], total, left) : 0;
      left -= out[i];
    }
  }
  return out;
}
