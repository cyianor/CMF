// Generated by cpp11: do not edit by hand
// clang-format off

#include <cpp11/R.hpp>
#include <Rcpp.h>
using namespace Rcpp;
#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// helper.cpp
void p_gradUsparse(const doubles_matrix<> Xm, doubles_matrix<> Gm_, const doubles_matrix<> CUm, const doubles_matrix<> OUm, const doubles_matrix<> Cm, const int idx, const double tau, const doubles Rowm, const doubles Colm);
extern "C" SEXP _CMF_p_gradUsparse(SEXP Xm, SEXP Gm_, SEXP CUm, SEXP OUm, SEXP Cm, SEXP idx, SEXP tau, SEXP Rowm, SEXP Colm) {
  BEGIN_CPP11
    p_gradUsparse(cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(Xm), cpp11::as_cpp<cpp11::decay_t<doubles_matrix<>>>(Gm_), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(CUm), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(OUm), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(Cm), cpp11::as_cpp<cpp11::decay_t<const int>>(idx), cpp11::as_cpp<cpp11::decay_t<const double>>(tau), cpp11::as_cpp<cpp11::decay_t<const doubles>>(Rowm), cpp11::as_cpp<cpp11::decay_t<const doubles>>(Colm));
    return R_NilValue;
  END_CPP11
}
// helper.cpp
void p_updatePseudoData(doubles_matrix<> Xm_, const doubles_matrix<> U1m, const doubles_matrix<> U2m, const doubles Rv, const doubles Cv);
extern "C" SEXP _CMF_p_updatePseudoData(SEXP Xm_, SEXP U1m, SEXP U2m, SEXP Rv, SEXP Cv) {
  BEGIN_CPP11
    p_updatePseudoData(cpp11::as_cpp<cpp11::decay_t<doubles_matrix<>>>(Xm_), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(U1m), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(U2m), cpp11::as_cpp<cpp11::decay_t<const doubles>>(Rv), cpp11::as_cpp<cpp11::decay_t<const doubles>>(Cv));
    return R_NilValue;
  END_CPP11
}
// helper.cpp
double p_updateTau(const doubles_matrix<> Xm, const doubles_matrix<> U1m, const doubles_matrix<> U2m, const doubles_matrix<> cov1m, const doubles_matrix<> cov2m, const doubles Rv, const doubles Cv, const doubles nu1v, const doubles nu2v);
extern "C" SEXP _CMF_p_updateTau(SEXP Xm, SEXP U1m, SEXP U2m, SEXP cov1m, SEXP cov2m, SEXP Rv, SEXP Cv, SEXP nu1v, SEXP nu2v) {
  BEGIN_CPP11
    return cpp11::as_sexp(p_updateTau(cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(Xm), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(U1m), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(U2m), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(cov1m), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(cov2m), cpp11::as_cpp<cpp11::decay_t<const doubles>>(Rv), cpp11::as_cpp<cpp11::decay_t<const doubles>>(Cv), cpp11::as_cpp<cpp11::decay_t<const doubles>>(nu1v), cpp11::as_cpp<cpp11::decay_t<const doubles>>(nu2v)));
  END_CPP11
}
// helper.cpp
sexp p_updateMean(const doubles_matrix<> Xm, const doubles_matrix<> U1m, const doubles_matrix<> U2m, const int idx, const doubles Mv);
extern "C" SEXP _CMF_p_updateMean(SEXP Xm, SEXP U1m, SEXP U2m, SEXP idx, SEXP Mv) {
  BEGIN_CPP11
    return cpp11::as_sexp(p_updateMean(cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(Xm), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(U1m), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(U2m), cpp11::as_cpp<cpp11::decay_t<const int>>(idx), cpp11::as_cpp<cpp11::decay_t<const doubles>>(Mv)));
  END_CPP11
}
// helper.cpp
void p_covUsparse(const doubles_matrix<> Xm, doubles_matrix<> Cm_, const doubles_matrix<> OUm, const doubles_matrix<> OCm, const int idx, const double tau);
extern "C" SEXP _CMF_p_covUsparse(SEXP Xm, SEXP Cm_, SEXP OUm, SEXP OCm, SEXP idx, SEXP tau) {
  BEGIN_CPP11
    p_covUsparse(cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(Xm), cpp11::as_cpp<cpp11::decay_t<doubles_matrix<>>>(Cm_), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(OUm), cpp11::as_cpp<cpp11::decay_t<const doubles_matrix<>>>(OCm), cpp11::as_cpp<cpp11::decay_t<const int>>(idx), cpp11::as_cpp<cpp11::decay_t<const double>>(tau));
    return R_NilValue;
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_CMF_p_covUsparse",       (DL_FUNC) &_CMF_p_covUsparse,       6},
    {"_CMF_p_gradUsparse",      (DL_FUNC) &_CMF_p_gradUsparse,      9},
    {"_CMF_p_updateMean",       (DL_FUNC) &_CMF_p_updateMean,       5},
    {"_CMF_p_updatePseudoData", (DL_FUNC) &_CMF_p_updatePseudoData, 5},
    {"_CMF_p_updateTau",        (DL_FUNC) &_CMF_p_updateTau,        9},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_CMF(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}