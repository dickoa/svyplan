#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP C_psu_menu_part(SEXP sizes, SEXP N, SEXP take, SEXP f_lo,
                            SEXP f_hi, SEXP cutoff, SEXP forced, SEXP m,
                            SEXP ord, SEXP f_stop, SEXP ld);
extern SEXP C_menu_group_min(SEXP x, SEXP ends);
extern SEXP C_menu_dominated(SEXP C, SEXP ord);
extern SEXP C_menu_step_down(SEXP cost, SEXP C, SEXP ends, SEXP pick,
                             SEXP ld);

static const R_CallMethodDef CallEntries[] = {
  {"C_psu_menu_part", (DL_FUNC) &C_psu_menu_part, 11},
  {"C_menu_group_min", (DL_FUNC) &C_menu_group_min, 2},
  {"C_menu_dominated", (DL_FUNC) &C_menu_dominated, 2},
  {"C_menu_step_down", (DL_FUNC) &C_menu_step_down, 5},
  {NULL, NULL, 0}
};

void attribute_visible R_init_svyplan(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
