#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <cpp11/matrix.hpp>
using namespace cpp11;

[[cpp11::register]] doubles_matrix<> p_gradUsparse(
    const integers_matrix<> indices, const doubles Xm, const doubles_matrix<> Gm,
    const doubles_matrix<> CUm, const doubles_matrix<> OUm, const doubles_matrix<> Cm,
    const int idx, const double tau, const doubles Rowm, const doubles Colm) {
  // double* const pGm = REAL(Gm.data());
  writable::doubles_matrix<> Gm_out(Gm);

  const int N = indices.nrow();
  const int K = Gm.ncol();

  // If idx == 1 => update a row entity
  // If idx == 2 => update a column entity
  if (idx == 1) {
    for (int n = 0; n < N; n++) {
      const int r = indices(n, 0) - 1;
      const int c = indices(n, 1) - 1;

      double tmp = 0.0;
      for (int k = 0; k < K; k++) {
        tmp += CUm(r, k) * OUm(c, k);
      }
      tmp += -Xm[n] + Rowm[r] + Colm[c];

      for (int k = 0; k < K; k++) {
        // pGm[r + k * K] += tau * (tmp * OUm(c, k) + CUm(r, k) * Cm(c, k));
        Gm_out(r, k) += tau * (tmp * OUm(c, k) + CUm(r, k) * Cm(c, k));
      }
    }
  } else {
    for (int n = 0; n < N; n++) {
      const int r = indices(n, 0) - 1;
      const int c = indices(n, 1) - 1;

      double tmp = 0.0;
      for (int k = 0; k < K; k++) {
        tmp += CUm(c, k) * OUm(r, k);
      }
      tmp += -Xm[n] + Rowm[r] + Colm[c];

      for (int k = 0; k < K; k++) {
        // pGm[c + k * K] += tau * (tmp * OUm(r, k) + CUm(c, k) * Cm(r, k));
        Gm_out(c, k) += tau * (tmp * OUm(r, k) + CUm(c, k) * Cm(r, k));
      }
    }
  }

  return Gm_out;
}

[[cpp11::register]] doubles p_updatePseudoData(const integers_matrix<> indices,
                                               const doubles_matrix<> U1m,
                                               const doubles_matrix<> U2m,
                                               const doubles Rv, const doubles Cv) {
  const int N = indices.nrow();
  const int K = U1m.ncol();

  writable::doubles out(N);

  for (int n = 0; n < N; n++) {
    const int r = indices(n, 0) - 1;
    const int c = indices(n, 1) - 1;

    double tmp = 0.0;
    for (int k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    out[n] = tmp + Rv[r] + Cv[c];
  }

  return out;
}

[[cpp11::register]] double p_updateTau(const integers_matrix<> indices, const doubles Xm,
                                       const doubles_matrix<> U1m,
                                       const doubles_matrix<> U2m,
                                       const doubles_matrix<> cov1m,
                                       const doubles_matrix<> cov2m, const doubles Rv,
                                       const doubles Cv, const doubles nu1v,
                                       const doubles nu2v) {
  const int N = indices.nrow();
  const int K = U1m.ncol();
  double out = 0.0;

  for (int n = 0; n < N; n++) {
    const int r = indices(n, 0) - 1;
    const int c = indices(n, 1) - 1;

    double tmp = 0.0;
    for (int k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    tmp += Rv[r] + Cv[c];
    tmp = Xm[n] - tmp;
    tmp = tmp * tmp;
    for (int k = 0; k < K; k++) {
      tmp += cov1m(r, k) * U2m(c, k) * U2m(c, k) + U1m(r, k) * U1m(r, k) * cov2m(c, k) +
             cov1m(r, k) * cov2m(c, k);
    }
    tmp += nu1v[r] + nu2v[c];

    out += tmp;
  }

  return out;
}

[[cpp11::register]] list p_updateMean(const integers_matrix<> indices, const doubles Xm,
                                      const doubles_matrix<> U1m,
                                      const doubles_matrix<> U2m, const int idx,
                                      const doubles Mv) {
  const int N = indices.nrow();
  const int K = U1m.ncol();

  // If idx == 1 => update a row entity
  // If idx == 2 => update a column entity
  writable::doubles Nv(idx == 1 ? U1m.nrow() : U2m.nrow());
  writable::integers Cv(idx == 1 ? U1m.nrow() : U2m.nrow());
  for (R_xlen_t i = 0; i < Cv.size(); i++) {
    Nv[i] = 0.0;
    Cv[i] = 0;
  }

  for (int n = 0; n < N; n++) {
    const int r = indices(n, 0) - 1;
    const int c = indices(n, 1) - 1;

    double tmp = 0.0;
    for (int k = 0; k < K; k++) {
      tmp += U1m(r, k) * U2m(c, k);
    }
    if (idx == 1) {
      tmp = Xm[n] - tmp - Mv[c];
      Nv[r] += tmp;
      Cv[r]++;
    } else {
      tmp = Xm[n] - tmp - Mv[r];
      Nv[c] += tmp;
      Cv[c]++;
    }
  }

  return writable::list({"sum"_nm = Nv, "count"_nm = Cv});
}

[[cpp11::register]] doubles_matrix<> p_covUsparse(const integers_matrix<> indices,
                                                  const doubles_matrix<> Cm,
                                                  const doubles_matrix<> OUm,
                                                  const doubles_matrix<> OCm,
                                                  const int idx, const double tau) {
  // double* const pCm = REAL(Cm.data());
  writable::doubles_matrix<> Cm_out(Cm);

  const int N = indices.nrow();
  const int K = Cm.ncol();

  if (idx == 1) {
    for (int n = 0; n < N; n++) {
      const int r = indices(n, 0) - 1;
      const int c = indices(n, 1) - 1;

      for (int k = 0; k < K; k++) {
        // pCm[r + k * K] += tau * (OUm(c, k) * OUm(c, k) + OCm(c, k));
        Cm_out(r, k) += tau * (OUm(c, k) * OUm(c, k) + OCm(c, k));
      }
    }
  } else {
    for (int n = 0; n < N; n++) {
      const int r = indices(n, 0) - 1;
      const int c = indices(n, 1) - 1;

      for (int k = 0; k < K; k++) {
        // pCm[c + k * K] += tau * (OUm(r, k) * OUm(r, k) + OCm(r, k));
        Cm_out(c, k) += tau * (OUm(r, k) * OUm(r, k) + OCm(r, k));
      }
    }
  }

  return Cm_out;
}
