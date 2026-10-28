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
  /* A vector gives vectors, which the node loop reads without a copy. */
  int vec = Rf_isNull(dim);
  SEXP mn = PROTECT(vec ? Rf_allocVector(REALSXP, G) :
                          Rf_allocMatrix(REALSXP, G, K));
  SEXP sec = PROTECT(vec ? Rf_allocVector(REALSXP, G) :
                           Rf_allocMatrix(REALSXP, G, K));
  SEXP arg = PROTECT(vec ? Rf_allocVector(INTSXP, G) :
                           Rf_allocMatrix(INTSXP, G, K));
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

/* Move each plan to cheaper states, largest saving first, without building
 * a residual-capacity matrix for every group. Groups are consecutive runs
 * ending at `ends`, and `pick` has one 1-based state per group. The
 * starting load is summed as .colSums() sums it, in long double when `ld`
 * is true, so every feasibility test matches the R scan's to the bit. */
SEXP C_menu_step_down(SEXP cost, SEXP C, SEXP ends, SEXP pick, SEXP ld) {
  if (!Rf_isReal(cost)) Rf_error("'cost' must be double");
  if (!Rf_isReal(C) || !Rf_isMatrix(C)) Rf_error("'C' must be a double matrix");
  if (!Rf_isInteger(ends) || !Rf_isInteger(pick)) {
    Rf_error("'ends' and 'pick' must be integer");
  }
  if (TYPEOF(ld) != LGLSXP || XLENGTH(ld) != 1) Rf_error("'ld' must be one logical");
  R_xlen_t rows = Rf_nrows(C), n = XLENGTH(cost);
  int K = Rf_ncols(C), H = LENGTH(ends);
  if (rows != n || LENGTH(pick) != H) Rf_error("menu dimensions do not agree");
  const int *e = INTEGER(ends), *ip = INTEGER(pick);
  for (int h = 0; h < H; h++) {
    R_xlen_t lo = h == 0 ? 0 : e[h - 1];
    if (e[h] <= lo || e[h] > rows) Rf_error("'ends' must increase within 'C'");
    if (ip[h] == NA_INTEGER || ip[h] < 1 || ip[h] > rows) {
      Rf_error("'pick' must index rows of 'C'");
    }
  }
  const double *pc = REAL(C), *pcost = REAL(cost);
  double *load = (double *) R_alloc((size_t) (K > 0 ? K : 1), sizeof(double));
  int use_ld = LOGICAL(ld)[0] == TRUE;
  for (int k = 0; k < K; k++) {
    const double *col = pc + (R_xlen_t) k * rows;
    if (use_ld) {
      long double sum = 0;
      for (int h = 0; h < H; h++) sum += col[ip[h] - 1];
      load[k] = (double) sum;
    } else {
      double sum = 0;
      for (int h = 0; h < H; h++) sum += col[ip[h] - 1];
      load[k] = sum;
    }
  }
  SEXP result = PROTECT(Rf_allocVector(INTSXP, H));
  int *out = INTEGER(result);
  for (int h = 0; h < H; h++) out[h] = ip[h];
  unsigned long scans = 0;
  for (;;) {
    int best_h = -1;
    R_xlen_t best_s = -1;
    double best_save = 0;
    for (int h = 0; h < H; h++) {
      if ((++scans & 1023UL) == 0) R_CheckUserInterrupt();
      R_xlen_t cur = out[h] - 1;
      R_xlen_t local = -1;
      for (R_xlen_t i = h == 0 ? 0 : e[h - 1]; i < e[h]; i++) {
        if (pcost[i] >= pcost[cur]) continue;
        int k = 0;
        while (k < K) {
          double residual = load[k] - pc[cur + (R_xlen_t) k * rows];
          if (pc[i + (R_xlen_t) k * rows] + residual > 1 + 1e-9) break;
          k++;
        }
        if (k == K && (local < 0 || pcost[i] < pcost[local])) local = i;
      }
      if (local >= 0) {
        double save = pcost[cur] - pcost[local];
        if (best_h < 0 || save > best_save) {
          best_h = h;
          best_s = local;
          best_save = save;
        }
      }
    }
    if (best_h < 0) break;
    R_xlen_t old = out[best_h] - 1;
    for (int k = 0; k < K; k++) {
      load[k] = load[k] - pc[old + (R_xlen_t) k * rows] +
        pc[best_s + (R_xlen_t) k * rows];
    }
    out[best_h] = (int) best_s + 1;
  }
  UNPROTECT(1);
  return result;
}
