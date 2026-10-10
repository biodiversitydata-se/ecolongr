#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// acc holds the pairs i < j in the order of a dist object: column i, rows j > i.
// Rows must be grouped by taxon; sample ids are 0-based.
//
// Within a taxon the samples are sorted, and each thread owns a range of sample
// ids. A pair is written by the thread that owns its smaller id, so no two
// threads write the same cell and every cell is added in taxon order whatever
// the number of threads.
// [[Rcpp::export]]
void accumulate_pairs(NumericVector acc, IntegerVector sample, IntegerVector taxon,
                      NumericVector quant, int n, bool product = false, int threads = 1) {
  const R_xlen_t N = sample.size();
  if (acc.size() != (R_xlen_t) n * (n - 1) / 2) stop("`acc` must have n * (n - 1) / 2 elements");
  for (R_xlen_t r = 0; r < N; ++r) {
    if (sample[r] < 0 || sample[r] >= n) stop("`sample` ids must be in [0, n)");
  }
  if (N == 0 || n < 2) return;
  if (threads < 1) threads = 1;

  std::vector<int> s(sample.begin(), sample.end());
  std::vector<double> q(quant.begin(), quant.end());
  std::vector<R_xlen_t> start;
  for (R_xlen_t r = 0; r < N; ++r) {
    if (r == 0 || taxon[r] != taxon[r - 1]) start.push_back(r);
  }
  start.push_back(N);
  const R_xlen_t G = (R_xlen_t) start.size() - 1;

  #pragma omp parallel for schedule(dynamic, 1024) num_threads(threads)
  for (R_xlen_t g = 0; g < G; ++g) {
    const R_xlen_t lo = start[g], hi = start[g + 1];
    if (std::is_sorted(s.begin() + lo, s.begin() + hi)) continue;
    std::vector<R_xlen_t> order(hi - lo);
    std::iota(order.begin(), order.end(), lo);
    std::sort(order.begin(), order.end(), [&](R_xlen_t x, R_xlen_t y) { return s[x] < s[y]; });
    std::vector<int> s2(hi - lo);
    std::vector<double> q2(hi - lo);
    for (R_xlen_t k = 0; k < hi - lo; ++k) { s2[k] = s[order[k]]; q2[k] = q[order[k]]; }
    std::copy(s2.begin(), s2.end(), s.begin() + lo);
    std::copy(q2.begin(), q2.end(), q.begin() + lo);
  }

  // Work per sample id: the number of larger ids sharing a taxon with it.
  std::vector<double> work(n, 0.0);
  for (R_xlen_t g = 0; g < G; ++g) {
    const R_xlen_t lo = start[g], hi = start[g + 1];
    for (R_xlen_t a = lo; a < hi; ++a) work[s[a]] += (double) (hi - 1 - a);
  }
  const int T = std::min(threads, n);
  std::vector<int> cut(T + 1, n);
  cut[0] = 0;
  const double total = std::accumulate(work.begin(), work.end(), 0.0);
  double seen = 0;
  int t = 1;
  for (int id = 0; id < n && t < T; ++id) {
    seen += work[id];
    while (t < T && seen >= total * t / T) cut[t++] = id + 1;
  }

  double* out = acc.begin();
  #pragma omp parallel for schedule(static, 1) num_threads(T)
  for (int part = 0; part < T; ++part) {
    const int idLo = cut[part], idHi = cut[part + 1];
    for (R_xlen_t g = 0; g < G; ++g) {
      const R_xlen_t lo = start[g], hi = start[g + 1];
      R_xlen_t a = std::lower_bound(s.begin() + lo, s.begin() + hi, idLo) - s.begin();
      for (; a < hi && s[a] < idHi; ++a) {
        const R_xlen_t i = s[a];
        double* row = out + i * (2 * (R_xlen_t) n - i - 1) / 2 - i - 1;
        const double qa = q[a];
        if (product) {
          for (R_xlen_t b = a + 1; b < hi; ++b) row[s[b]] += qa * q[b];
        } else {
          for (R_xlen_t b = a + 1; b < hi; ++b) row[s[b]] += qa < q[b] ? qa : q[b];
        }
      }
    }
  }
}

// method: 0 = bray, 1 = jaccard, 2 = euclidean. acc holds the shared minima
// (or products for euclidean) and is replaced by the distances.
// [[Rcpp::export]]
void finalize_dist(NumericVector acc, NumericVector total, NumericVector square,
                   int n, int method, int threads = 1) {
  if (acc.size() != (R_xlen_t) n * (n - 1) / 2) stop("`acc` must have n * (n - 1) / 2 elements");
  double* out = acc.begin();
  #pragma omp parallel for schedule(dynamic, 16) num_threads(threads)
  for (int i = 0; i < n - 1; ++i) {
    R_xlen_t k = (R_xlen_t) i * (2 * (R_xlen_t) n - i - 1) / 2;
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
