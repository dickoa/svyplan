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
extern SEXP C_stage_search(SEXP R_sexp, SEXP rates_sexp, SEXP M_sexp,
                           SEXP cv2_target_sexp, SEXP cost_sexp,
                           SEXP lower_sexp, SEXP upper_sexp,
                           SEXP budget_var_sexp, SEXP min_total_sexp,
                           SEXP n_init_sexp, SEXP max_nodes_sexp);
extern SEXP C_menu_tree(SEXP cost_s, SEXP C_s, SEXP ends_s, SEXP lambda_s,
                        SEXP start_s, SEXP grain_s, SEXP max_nodes_s,
                        SEXP max_open_s, SEXP gap_s, SEXP exact_s);

static const R_CallMethodDef CallEntries[] = {
  {"C_psu_menu_part", (DL_FUNC) &C_psu_menu_part, 11},
  {"C_menu_group_min", (DL_FUNC) &C_menu_group_min, 2},
  {"C_menu_dominated", (DL_FUNC) &C_menu_dominated, 2},
  {"C_menu_step_down", (DL_FUNC) &C_menu_step_down, 5},
  {"C_stage_search", (DL_FUNC) &C_stage_search, 11},
  {"C_menu_tree", (DL_FUNC) &C_menu_tree, 10},
  {NULL, NULL, 0}
};

void attribute_visible R_init_svyplan(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
