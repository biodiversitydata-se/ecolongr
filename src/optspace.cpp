#include <Rcpp.h>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;

// Observed cells are given as (i, j, value) triplets with 0-based indices. The factor
// matrices have r columns; every function loops over the cells once.

// Fitted value of every cell: the dot product of row i of XS and row j of Y.
// [[Rcpp::export]]
NumericVector cell_fitted(NumericMatrix XS, NumericMatrix Y, IntegerVector i, IntegerVector j, int threads = 1) {
  const R_xlen_t N = i.size();
  const int r = XS.ncol();
  const int n = XS.nrow(), p = Y.nrow();
  NumericVector out(N);
  #pragma omp parallel for num_threads(threads) schedule(static)
  for (R_xlen_t c = 0; c < N; ++c) {
    double s = 0;
    for (int a = 0; a < r; ++a) s += XS[i[c] + (R_xlen_t) a * n] * Y[j[c] + (R_xlen_t) a * p];
    out[c] = s;
  }
  return out;
}

// out[row_c, ] += val_c * M[col_c, ] for every cell: a sparse matrix times a dense one.
// [[Rcpp::export]]
NumericMatrix cell_multiply(IntegerVector row, IntegerVector col, NumericVector val, NumericMatrix M,
                            int nrow_out, int threads = 1) {
  const R_xlen_t N = row.size();
  const int r = M.ncol();
  const int m = M.nrow();
  const R_xlen_t size = (R_xlen_t) nrow_out * r;
  if (threads < 1) threads = 1;
  std::vector<std::vector<double> > local(threads, std::vector<double>(size, 0.0));
  #pragma omp parallel num_threads(threads)
  {
    int t = 0;
    #ifdef _OPENMP
    t = omp_get_thread_num();
    #endif
    double* acc = local[t].data();
    #pragma omp for schedule(static)
    for (R_xlen_t c = 0; c < N; ++c) {
      const double v = val[c];
      for (int a = 0; a < r; ++a) acc[row[c] + (R_xlen_t) a * nrow_out] += v * M[col[c] + (R_xlen_t) a * m];
    }
  }
  NumericMatrix out(nrow_out, r);
  for (int t = 0; t < threads; ++t) for (R_xlen_t k = 0; k < size; ++k) out[k] += local[t][k];
  return out;
}

// The r^2 x r^2 matrix A with A[b * r + a, j * r + i] = sum over cells of X[k, a] X[k, i] Y[l, j] Y[l, b].
// [[Rcpp::export]]
NumericMatrix cell_gram(NumericMatrix X, NumericMatrix Y, IntegerVector i, IntegerVector j, int threads = 1) {
  const R_xlen_t N = i.size();
  const int r = X.ncol();
  const int n = X.nrow(), p = Y.nrow();
  const int rr = r * r;
  if (threads < 1) threads = 1;
  std::vector<std::vector<double> > local(threads, std::vector<double>((size_t) rr * rr, 0.0));
  #pragma omp parallel num_threads(threads)
  {
    int t = 0;
    #ifdef _OPENMP
    t = omp_get_thread_num();
    #endif
    double* acc = local[t].data();
    std::vector<double> xk(r), yl(r);
    #pragma omp for schedule(static)
    for (R_xlen_t c = 0; c < N; ++c) {
      for (int a = 0; a < r; ++a) { xk[a] = X[i[c] + (R_xlen_t) a * n]; yl[a] = Y[j[c] + (R_xlen_t) a * p]; }
      for (int jj = 0; jj < r; ++jj) for (int ii = 0; ii < r; ++ii) {
        const int col = jj * r + ii;
        for (int b = 0; b < r; ++b) for (int a = 0; a < r; ++a) {
          acc[(size_t) col * rr + (b * r + a)] += xk[a] * xk[ii] * yl[jj] * yl[b];
        }
      }
    }
  }
  NumericMatrix out(rr, rr);
  for (int t = 0; t < threads; ++t) for (size_t k = 0; k < (size_t) rr * rr; ++k) out[k] += local[t][k];
  return out;
}
