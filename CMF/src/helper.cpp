#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void p_gradUsparse(NumericMatrix Xm, NumericMatrix Gm, NumericMatrix CUm,
                   NumericMatrix OUm, NumericMatrix Cm, int idx,
                   double tau, NumericVector Rowm, NumericVector Colm) {
  const size_t N = static_cast<size_t>(Xm.nrow());
  const size_t K = static_cast<size_t>(Gm.ncol());

  // If idx == 1 => update a row entity
  // If idx == 2 => update a column entity
  if (idx == 1) {
    for (size_t n = 0; n < N; n++) {
      const size_t r = static_cast<size_t>(Xm(n, 0) - 1);
      const size_t c = static_cast<size_t>(Xm(n, 1) - 1);
      double tmp = 0.0;
      for (size_t k = 0; k < K; k++) {
        tmp += CUm(r, k) * OUm(c, k);
      }
      tmp += -Xm(n, 2) + Rowm(r) + Colm(c);
      for (int k = 0; k < K; k++) {
        Gm(r, k) += tau * (tmp * OUm(c, k) + CUm(r, k) * Cm(c, k));
      }
    }
  }
  else {
    for (size_t n = 0; n < N; n++) {
      const size_t r = static_cast<size_t>(Xm(n, 0) - 1);
      const size_t c = static_cast<size_t>(Xm(n, 1) - 1);
      double tmp = 0.0;
      for (size_t k = 0; k < K; k++) {
        tmp += CUm(c, k) * OUm(r, k);
      }
      tmp += -Xm(n, 2) + Rowm(r) + Colm(c);
      for (size_t k = 0; k < K; k++) {
        Gm(c, k) += tau * (tmp * OUm(r, k) + CUm(c, k) * Cm(r, k));
      }
    }
  }
}

// [[Rcpp::export]]
void p_updatePseudoData(NumericMatrix Xm, NumericMatrix U1m,
                        NumericMatrix U2m, NumericVector Rv,
                        NumericVector Cv) {
  const size_t N = static_cast<size_t>(Xm.nrow());
  const size_t K = static_cast<size_t>(U1m.ncol());

  for (size_t n = 0; n < N; n++) {
    const size_t r = static_cast<size_t>(Xm(n, 0) - 1);
    const size_t c = static_cast<size_t>(Xm(n, 1) - 1);
    double tmp = 0.0;
    for (size_t k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    tmp += Rv(r) + Cv(c);
    Xm(n, 2) = tmp;
  }
}

// [[Rcpp::export]]
double p_updateTau(NumericMatrix Xm, NumericMatrix U1m,
                   NumericMatrix U2m, NumericMatrix cov1m,
                   NumericMatrix cov2m, NumericVector Rv,
                   NumericVector Cv, NumericVector nu1v,
                   NumericVector nu2v) {
  const size_t N = static_cast<size_t>(Xm.nrow());
  const size_t K = static_cast<size_t>(U1m.ncol());
  double out = 0.0;

  for (size_t n = 0; n < N; n++) {
    const size_t r = static_cast<size_t>(Xm(n, 0) - 1);
    const size_t c = static_cast<size_t>(Xm(n, 1) - 1);
    double tmp = 0.0;
    for (size_t k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    tmp += Rv(r) + Cv(c);
    tmp = Xm(n, 2) - tmp;
    tmp = tmp * tmp;
    for (size_t k = 0; k < K; k++) {
      tmp += cov1m(r, k) * U2m(c, k) * U2m(c, k)
             + U1m(r, k) * U1m(r, k) * cov2m(c, k)
             + cov1m(r, k) * cov2m(c, k);
    }
    tmp += nu1v(r) + nu2v(c);
    out += tmp;
  }

  return(out);
}

// [[Rcpp::export]]
List p_updateMean(NumericMatrix Xm, NumericMatrix U1m, NumericMatrix U2m,
                  int idx, NumericVector Mv) {
  const size_t N = static_cast<size_t>(Xm.nrow());
  const size_t K = static_cast<size_t>(U1m.ncol());
  NumericVector Nv, Cv;

  // If idx == 1 => update a row entity
  // If idx == 2 => update a column entity
  if (idx == 1) {
    Nv = NumericVector(U1m.nrow());
    Cv = NumericVector(U1m.nrow());
  } else {
    Nv = NumericVector(U2m.nrow());
    Cv = NumericVector(U2m.nrow());
  }

  for (size_t n = 0; n < N; n++) {
    const size_t r = static_cast<size_t>(Xm(n, 0) - 1);
    const size_t c = static_cast<size_t>(Xm(n, 1) - 1);
    double tmp = 0.0;
    for (size_t k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    if (idx == 1) {
      tmp = Xm(n, 2) - tmp - Mv(c);
      Nv(r) += tmp;
      Cv(r)++;
    } else {
      tmp = Xm(n, 2) - tmp - Mv(r);
      Nv(c) += tmp;
      Cv(c)++;
    }
  }

  List ret = List::create(_["sum"] = Nv, _["count"] = Cv);
  return(ret);
}

// [[Rcpp::export]]
void p_covUsparse(NumericMatrix Xm, NumericMatrix Cm, NumericMatrix OUm,
                  NumericMatrix OCm, int idx, double tau) {
  const size_t N = static_cast<size_t>(Xm.nrow());
  const size_t K = static_cast<size_t>(Cm.ncol());

  if (idx == 1) {
    for (size_t n = 0; n < N; n++) {
      const size_t r = static_cast<size_t>(Xm(n, 0) - 1);
      const size_t c = static_cast<size_t>(Xm(n, 1) - 1);
      for (size_t k = 0; k < K; k++) {
        Cm(r, k) += tau * (OUm(c, k) * OUm(c, k) + OCm(c, k));
      }
    }
  } else {
    for (size_t n = 0; n < N; n++) {
      const size_t r = static_cast<size_t>(Xm(n, 0) - 1);
      const size_t c = static_cast<size_t>(Xm(n, 1) - 1);
      for (size_t k = 0; k < K; k++) {
        Cm(c, k) += tau * (OUm(r, k) * OUm(r, k) + OCm(r, k));
      }
    }
  }
}
