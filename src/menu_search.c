/* Group minima for the menu search.
 *
 * The search keeps a node's states in consecutive runs, one per group, and
 * needs for each run and column the smallest value, its position and the
 * second smallest. R has no single call for that, and a loop over groups
 * dominates the time of a node. Ties go to the first position, as in
 * which.min(). */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

SEXP C_menu_group_min(SEXP x, SEXP ends) {
  if (!Rf_isReal(x)) Rf_error("'x' must be double");
  if (!Rf_isInteger(ends)) Rf_error("'ends' must be integer");
  R_xlen_t n = XLENGTH(x);
  int G = LENGTH(ends);
  SEXP dim = Rf_getAttrib(x, R_DimSymbol);
  R_xlen_t rows = n;
  int K = 1;
  if (!Rf_isNull(dim)) {
    if (LENGTH(dim) != 2) Rf_error("'x' must be a vector or a matrix");
    rows = INTEGER(dim)[0];
    K = INTEGER(dim)[1];
  }
  const int *e = INTEGER(ends);
  for (int g = 0; g < G; g++) {
    R_xlen_t lo = g == 0 ? 0 : e[g - 1];
    if (e[g] <= lo || e[g] > rows) Rf_error("'ends' must increase within 'x'");
  }
  SEXP mn = PROTECT(Rf_allocMatrix(REALSXP, G, K));
  SEXP sec = PROTECT(Rf_allocMatrix(REALSXP, G, K));
  SEXP arg = PROTECT(Rf_allocMatrix(INTSXP, G, K));
  const double *px = REAL(x);
  double *pm = REAL(mn), *ps = REAL(sec);
  int *pa = INTEGER(arg);
  for (int k = 0; k < K; k++) {
    const double *col = px + (R_xlen_t) k * rows;
    for (int g = 0; g < G; g++) {
      R_xlen_t lo = g == 0 ? 0 : e[g - 1];
      double m1 = col[lo], m2 = R_PosInf;
      R_xlen_t a = lo;
      for (R_xlen_t i = lo + 1; i < e[g]; i++) {
        double v = col[i];
        if (v < m1) {
          m2 = m1;
          m1 = v;
          a = i;
        } else if (v < m2) {
          m2 = v;
        }
      }
      R_xlen_t at = g + (R_xlen_t) k * G;
      pm[at] = m1;
      ps[at] = m2;
      pa[at] = (int) (a + 1);
    }
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, mn);
  SET_VECTOR_ELT(out, 1, arg);
  SET_VECTOR_ELT(out, 2, sec);
  SEXP nm = PROTECT(Rf_allocVector(STRSXP, 3));
  SET_STRING_ELT(nm, 0, Rf_mkChar("min"));
  SET_STRING_ELT(nm, 1, Rf_mkChar("arg"));
  SET_STRING_ELT(nm, 2, Rf_mkChar("second"));
  Rf_setAttrib(out, R_NamesSymbol, nm);
  UNPROTECT(5);
  return out;
}

/* Dominance among one group's states.
 *
 * `ord` lists states (rows of `C`, from 1) in the order they are visited,
 * cost first. A state is dropped when a state kept before it loads no
 * column more, read as C[kept, k] <= C[i, k] + 1e-15 in every column, as
 * the R loop it replaces compared them. The R loop held every kept state
 * against every later one and took most of a root on registers with
 * thousands of states per stratum. Returns, for each entry of `ord`,
 * whether it is kept. */
SEXP C_menu_dominated(SEXP C, SEXP ord) {
  if (!Rf_isReal(C) || !Rf_isMatrix(C)) Rf_error("'C' must be a double matrix");
  if (!Rf_isInteger(ord)) Rf_error("'ord' must be integer");
  R_xlen_t rows = Rf_nrows(C);
  int K = Rf_ncols(C);
  int m = LENGTH(ord);
  const int *o = INTEGER(ord);
  for (int a = 0; a < m; a++) {
    if (o[a] == NA_INTEGER || o[a] < 1 || o[a] > rows) {
      Rf_error("'ord' must index rows of 'C'");
    }
  }
  const double *pc = REAL(C);
  SEXP out = PROTECT(Rf_allocVector(LGLSXP, m));
  int *keep = LOGICAL(out);
  R_xlen_t *kept = (R_xlen_t *) R_alloc((size_t) (m > 0 ? m : 1),
                                        sizeof(R_xlen_t));
  int nk = 0;
  for (int a = 0; a < m; a++) {
    if ((a & 1023) == 0) R_CheckUserInterrupt();
    R_xlen_t i = o[a] - 1;
    int dominated = 0;
    for (int b = 0; b < nk && !dominated; b++) {
      R_xlen_t j = kept[b];
      int k = 0;
      while (k < K && pc[j + k * rows] <= pc[i + k * rows] + 1e-15) k++;
      dominated = k == K;
    }
    keep[a] = !dominated;
    if (!dominated) kept[nk++] = i;
  }
  UNPROTECT(1);
  return out;
}
