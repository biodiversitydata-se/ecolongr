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
// Within a taxon the samples are sorted, and each block of sample ids is handled
// by one thread at a time. A pair is written in the block that owns its smaller
// id, so no two threads write the same cell and every cell is added in taxon
// order whatever the number of threads.
// [[Rcpp::export]]
void accumulate_pairs(NumericVector acc, IntegerVector sample, IntegerVector taxon,
                      NumericVector quant, int n, bool product = false, int threads = 1,
                      bool canberra = false, NumericVector count = NumericVector(0)) {
  const R_xlen_t N = sample.size();
  if (acc.size() != (R_xlen_t) n * (n - 1) / 2) stop("`acc` must have n * (n - 1) / 2 elements");
  for (R_xlen_t r = 0; r < N; ++r) {
    if (sample[r] < 0 || sample[r] >= n) stop("`sample` ids must be in [0, n)");
  }
  const bool counting = canberra && count.size() > 0;
  if (counting && count.size() != acc.size()) stop("`count` must have the length of `acc`");
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

  // Blocks of consecutive sample ids, each owning about 1 MB of acc so that it stays
  // in cache while all taxa are applied to it; at least four per thread for balance.
  const R_xlen_t cells = (R_xlen_t) n * (n - 1) / 2;
  const R_xlen_t blockCells = 131072;
  int P = (int) std::min<R_xlen_t>(n, std::max<R_xlen_t>(4 * (R_xlen_t) threads, (cells + blockCells - 1) / blockCells));
  std::vector<int> cut(P + 1, n);
  cut[0] = 0;
  {
    int part = 1;
    R_xlen_t seen = 0;
    for (int id = 0; id < n && part < P; ++id) {
      seen += n - 1 - id;
      while (part < P && seen >= cells * part / P) cut[part++] = id + 1;
    }
  }

  R_xlen_t maxK = 0;
  for (R_xlen_t g = 0; g < G; ++g) maxK = std::max(maxK, start[g + 1] - start[g]);

  double* out = acc.begin();
  double* counts = counting ? count.begin() : nullptr;
  #pragma omp parallel num_threads(threads)
  {
    std::vector<double> buffer(maxK);
    double* __restrict tmp = buffer.data();
    #pragma omp for schedule(dynamic, 1)
    for (int part = 0; part < P; ++part) {
      const int idLo = cut[part], idHi = cut[part + 1];
      if (idLo >= idHi) continue;
      for (R_xlen_t g = 0; g < G; ++g) {
        const R_xlen_t lo = start[g], hi = start[g + 1];
        R_xlen_t a = std::lower_bound(s.begin() + lo, s.begin() + hi, idLo) - s.begin();
        for (; a < hi && s[a] < idHi; ++a) {
          const R_xlen_t i = s[a];
          const R_xlen_t offset = i * (2 * (R_xlen_t) n - i - 1) / 2 - i - 1;
          double* __restrict row = out + offset;
          const double qa = q[a];
          const int* __restrict sb = s.data() + a + 1;
          const double* __restrict qb = q.data() + a + 1;
          const R_xlen_t len = hi - a - 1;
          if (canberra) {
            for (R_xlen_t j = 0; j < len; ++j) tmp[j] = std::fabs(qa - qb[j]) / (qa + qb[j]);
          } else if (product) {
            for (R_xlen_t j = 0; j < len; ++j) tmp[j] = qa * qb[j];
          } else {
            for (R_xlen_t j = 0; j < len; ++j) tmp[j] = qa < qb[j] ? qa : qb[j];
          }
          for (R_xlen_t j = 0; j < len; ++j) row[sb[j]] += tmp[j];
          if (counting) {
            double* __restrict crow = counts + offset;
            for (R_xlen_t j = 0; j < len; ++j) crow[sb[j]] += 1.0;
          }
        }
      }
    }
  }
}

// method: 0 = bray, 1 = jaccard, 2 = euclidean, 3 = aitchison, 4 = manhattan, 5 = canberra.
// acc holds the shared minima (products for euclidean and aitchison, the sum of
// |a - b| / (a + b) over shared taxa for canberra) and is replaced by the distances.
// For canberra, nnz holds the number of non-zero cells per sample and shared the number of
// taxa each pair has in common.
// For aitchison, total holds the row sums of the non-constant part of the clr and ntaxa
// the number of taxa: the clr is that part minus its row mean.
// [[Rcpp::export]]
void finalize_dist(NumericVector acc, NumericVector total, NumericVector square,
                   int n, int method, int threads = 1, double ntaxa = 1,
                   NumericVector nnz = NumericVector(0), NumericVector shared = NumericVector(0)) {
  if (acc.size() != (R_xlen_t) n * (n - 1) / 2) stop("`acc` must have n * (n - 1) / 2 elements");
  double* out = acc.begin();
  #pragma omp parallel for schedule(dynamic, 16) num_threads(threads)
  for (int i = 0; i < n - 1; ++i) {
    R_xlen_t k = (R_xlen_t) i * (2 * (R_xlen_t) n - i - 1) / 2;
    for (int j = i + 1; j < n; ++j, ++k) {
      if (method == 4) {
        const double d = total[i] + total[j] - 2 * out[k];
        out[k] = d > 0 ? d : 0.0;
      } else if (method == 5) {
        const double both = shared[k];
        const double union_size = nnz[i] + nnz[j] - both;
        out[k] = union_size > 0 ? ((nnz[i] - both) + (nnz[j] - both) + out[k]) / union_size : NA_REAL;
      } else if (method == 2 || method == 3) {
        double d2 = square[i] + square[j] - 2 * out[k];
        if (method == 3) d2 -= (total[i] - total[j]) * (total[i] - total[j]) / ntaxa;
        out[k] = d2 > 0 ? std::sqrt(d2) : 0.0;
      } else {
        const double b = 1 - 2 * out[k] / (total[i] + total[j]);
        out[k] = method == 1 ? 2 * b / (1 + b) : b;
      }
    }
  }
}
