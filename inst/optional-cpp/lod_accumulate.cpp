#include <Rcpp.h>
using namespace Rcpp;
// Row-major accumulation: inner loop writes contiguously.
// [[Rcpp::export]]
NumericMatrix lod_acc2(IntegerMatrix ocode, IntegerMatrix ccode,
                       NumericVector LT, int G) {
  int no = ocode.nrow(), nc = ccode.nrow(), nl = ocode.ncol();
  const double* lt = REAL(LT);
  std::vector<int> cb((size_t)nl * nc);
  for (int l = 0; l < nl; ++l)
    for (int j = 0; j < nc; ++j) {
      int c = ccode(j, l); if (c == NA_INTEGER) c = G;
      cb[(size_t)l * nc + j] = (c - 1) * G + l * G * G;
    }
  std::vector<double> buf((size_t)nc);
  NumericMatrix out(no, nc);
  for (int i = 0; i < no; ++i) {
    std::fill(buf.begin(), buf.end(), 0.0);
    for (int l = 0; l < nl; ++l) {
      int a = ocode(i, l); if (a == NA_INTEGER) a = G;
      const int off = a - 1;
      const int* cbl = &cb[(size_t)l * nc];
      double* b = buf.data();
      for (int j = 0; j < nc; ++j) b[j] += lt[cbl[j] + off];
    }
    for (int j = 0; j < nc; ++j) out(i, j) = buf[j];
  }
  return out;
}
