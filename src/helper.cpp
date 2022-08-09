#include <cpp11/matrix.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/list.hpp>
using namespace cpp11;

[[cpp11::register]] void p_gradUsparse(const doubles_matrix<> Xm,
                                       doubles_matrix<> Gm_,
                                       const doubles_matrix<> CUm,
                                       const doubles_matrix<> OUm,
                                       const doubles_matrix<> Cm,
                                       const int idx,
                                       const double tau,
                                       const doubles Rowm,
                                       const doubles Colm) {
  writable::doubles_matrix<> Gm(std::move(Gm_));

  const int N = Xm.nrow();
  const int K = Gm.ncol();

  // If idx == 1 => update a row entity
  // If idx == 2 => update a column entity
  if (idx == 1) {
    for (int n = 0; n < N; n++) {
      const int r = Xm(n, 0) - 1;
      const int c = Xm(n, 1) - 1;
      double tmp = 0.0;

      for (int k = 0; k < K; k++) {
        tmp += CUm(r, k) * OUm(c, k);
      }
      tmp += -Xm(n, 2) + Rowm[r] + Colm[c];

      for (int k = 0; k < K; k++) {
        Gm(r, k) += tau * (tmp * OUm(c, k) + CUm(r, k) * Cm(c, k));
      }
    }
  } else {
    for (int n = 0; n < N; n++) {
      const int r = Xm(n, 0) - 1;
      const int c = Xm(n, 1) - 1;
      double tmp = 0.0;

      for (int k = 0; k < K; k++) {
        tmp += CUm(c, k) * OUm(r, k);
      }
      tmp += -Xm(n, 2) + Rowm[r] + Colm[c];

      for (int k = 0; k < K; k++) {
        Gm(c, k) += tau * (tmp * OUm(r, k) + CUm(c, k) * Cm(r, k));
      }
    }
  }
}

[[cpp11::register]] void p_updatePseudoData(doubles_matrix<> Xm_,
                                            const doubles_matrix<> U1m,
                                            const doubles_matrix<> U2m,
                                            const doubles Rv,
                                            const doubles Cv) {
  writable::doubles_matrix<> Xm(std::move(Xm_));

  const int N = Xm.nrow();
  const int K = U1m.ncol();

  for (int n = 0; n < N; n++) {
    const int r = Xm(n, 0) - 1;
    const int c = Xm(n, 1) - 1;
    double tmp = 0.0;

    for (int k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    tmp += Rv[r] + Cv[c];

    Xm(n, 2) = tmp;
  }
}

[[cpp11::register]] double p_updateTau(const doubles_matrix<> Xm,
                                       const doubles_matrix<> U1m,
                                       const doubles_matrix<> U2m,
                                       const doubles_matrix<> cov1m,
                                       const doubles_matrix<> cov2m,
                                       const doubles Rv,
                                       const doubles Cv,
                                       const doubles nu1v,
                                       const doubles nu2v) {
  const int N = Xm.nrow();
  const int K = U1m.ncol();
  double out = 0.0;

  for (int n = 0; n < N; n++) {
    const int r = Xm(n, 0) - 1;
    const int c = Xm(n, 1) - 1;
    double tmp = 0.0;

    for (int k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    tmp += Rv[r] + Cv[c];
    tmp = Xm(n, 2) - tmp;
    tmp = tmp * tmp;
    for (int k = 0; k < K; k++) {
      tmp += cov1m(r, k) * U2m(c, k) * U2m(c, k)
             + U1m(r, k) * U1m(r, k) * cov2m(c, k)
             + cov1m(r, k) * cov2m(c, k);
    }
    tmp += nu1v[r] + nu2v[c];

    out += tmp;
  }

  return out;
}

[[cpp11::register]] sexp p_updateMean(const doubles_matrix<> Xm,
                                      const doubles_matrix<> U1m,
                                      const doubles_matrix<> U2m,
                                      const int idx,
                                      const doubles Mv) {
  const int N = Xm.nrow();
  const int K = U1m.ncol();

  // If idx == 1 => update a row entity
  // If idx == 2 => update a column entity
  writable::doubles Nv(idx == 1 ? U1m.nrow() : U2m.nrow());
  writable::integers Cv(idx == 1 ? U1m.nrow() : U2m.nrow());
  for (int i = 0; i < Cv.size(); i++) {
    Nv[i] = 0.0;
    Cv[i] = 0;
  }

  for (int n = 0; n < N; n++) {
    const int r = Xm(n, 0) - 1;
    const int c = Xm(n, 1) - 1;
    double tmp = 0.0;

    for (int k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }

    if (idx == 1) {
      tmp = Xm(n, 2) - tmp - Mv[c];
      Nv[r] += tmp;
      Cv[r]++;
    } else {
      tmp = Xm(n, 2) - tmp - Mv[r];
      Nv[c] += tmp;
      Cv[c]++;
    }
  }

  return writable::list({
    "sum"_nm = Nv, "count"_nm = Cv
  });
}

[[cpp11::register]] void p_covUsparse(const doubles_matrix<> Xm,
                                      doubles_matrix<> Cm_,
                                      const doubles_matrix<> OUm,
                                      const doubles_matrix<> OCm,
                                      const int idx,
                                      const double tau) {
  writable::doubles_matrix<> Cm(std::move(Cm_));

  const int N = Xm.nrow();
  const int K = Cm.ncol();

  if (idx == 1) {
    for (int n = 0; n < N; n++) {
      const int r = Xm(n, 0) - 1;
      const int c = Xm(n, 1) - 1;

      for (int k = 0; k < K; k++) {
        Cm(r, k) += tau * (OUm(c, k) * OUm(c, k) + OCm(c, k));
      }
    }
  } else {
    for (int n = 0; n < N; n++) {
      const int r = Xm(n, 0) - 1;
      const int c = Xm(n, 1) - 1;

      for (int k = 0; k < K; k++) {
        Cm(c, k) += tau * (OUm(r, k) * OUm(r, k) + OCm(r, k));
      }
    }
  }
}
