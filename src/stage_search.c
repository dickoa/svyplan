#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <float.h>

/*
 * Exact branch-and-bound whole-unit solver for S-stage cluster designs (S >= 2)
 * with J >= 1 indicators, optional stage-by-stage finite population correction,
 * stage-specific response rates, box constraints (covering fixed stages), and
 * optional minimum total sample size.
 *
 * All arrays are 0-indexed. Matrices (R, rates, M) are J x S in R column-major
 * order: element (j, s) is at index j + s * J.
 */

typedef struct {
  int S;
  int J;
  int budget_mode;
  double budget_var;
  double min_total;
  double max_nodes;
  double nodes;
  int hit_limit;

  const double *R;          /* J x S */
  const double *rates;      /* J x S */
  const double *M;          /* J x S (Inf = no FPC at stage s) */
  const double *cv2_target; /* J */
  const double *cost;       /* S */
  const double *lower;      /* S (whole numbers >= 1) */
  const double *upper;      /* S (whole numbers or Inf) */
  const double *n_init;     /* S (continuous warm start, or NULL) */

  double *rbar;             /* J x S: cumulative product of rates */
  double *alpha;            /* J x S: telescoping coefficients */
  double *kappa;            /* J: R[j,0] / M[j,0] (0 if M[j,0] is Inf) */
  double *W;                /* J x S: Cauchy-Schwarz suffix weight for stages > s */
  int    *suffix_pos;       /* J x S: 1 if all alpha[j, k] >= 0 for k > s */
  double *g_tail;           /* J x S: min variance tail at upper[s+1..S-1] */
  double *rem_cost_min;     /* S: min cost per unit of P_s for stages s..S-1 */

  /* Current recursion path */
  double *cur_n;            /* S */
  double *cur_P;            /* S */
  double *cur_C;            /* S */
  double *cur_fpc_sum;      /* J x S: partial sum of positive stage-FPC terms 0..s */
  double *cur_tel_sum;      /* J x S: partial sum of telescoping terms 0..s minus kappa */

  /* Best feasible incumbent */
  int found;
  long n_improve;           /* incumbent updates, read by the warm-start polish */
  double best_prim;         /* max ratio (budget mode) or cost (cv mode) */
  double best_tn;           /* total sample size, cv-mode tie-break */
  double best_take;         /* sum of takes, last tie-break */
  double best_cost;
  double best_ratio;
  double *best_n;           /* S */
} StageCtx;

static const double STAGE_TOL = 1e-9;
static const double INT_MAX_D = 2147483647.0;

/* Tolerances are relative to the quantity they guard, so a problem gives
 * the same design whatever its cost or variance scale. */
static inline double rel_slack(double ref) {
  return 1e-9 * fmax(fabs(ref), DBL_MIN);
}

/* The budget carries rounding only. A relative slack would buy whole units
 * once the budget is a billion times a unit cost. */
static inline double budget_cap(const StageCtx *ctx) {
  return ctx->budget_var + 8.0 * DBL_EPSILON * fabs(ctx->budget_var);
}

static inline double cost_cap(const StageCtx *ctx) {
  return ctx->best_cost + rel_slack(ctx->best_cost);
}

static inline double ratio_cap(const StageCtx *ctx) {
  return ctx->best_ratio + rel_slack(ctx->best_ratio);
}

static inline int cv_over(double cv2, double target) {
  return cv2 > target * (1.0 + 1e-9);
}

/* Lexicographic order: the primary objective within a relative tie
 * tolerance, then (cv mode) the total size, then the sum of takes. */
static int better(const StageCtx *ctx, double prim, double tn, double ts) {
  if (!ctx->found) return 1;
  double tol = 1e-12 * fmax(fmax(fabs(prim), fabs(ctx->best_prim)), DBL_MIN);
  if (prim < ctx->best_prim - tol) return 1;
  if (prim > ctx->best_prim + tol) return 0;
  if (!ctx->budget_mode) {
    if (tn < ctx->best_tn - 1e-9) return 1;
    if (tn > ctx->best_tn + 1e-9) return 0;
  }
  return ts < ctx->best_take - 1e-9;
}

static void store_incumbent(StageCtx *ctx, double prim, double tn, double ts,
                            double cost, double ratio, double n_last) {
  int last = ctx->S - 1;
  ctx->found = 1;
  ctx->n_improve++;
  ctx->best_prim = prim;
  ctx->best_tn = tn;
  ctx->best_take = ts;
  ctx->best_cost = cost;
  ctx->best_ratio = ratio;
  for (int s = 0; s < last; s++) ctx->best_n[s] = ctx->cur_n[s];
  ctx->best_n[last] = n_last;
}

/* Exact check of a whole last-stage size, in the positive stage-FPC form */
static int leaf_cv_ok(const StageCtx *ctx, int last, double P_prev,
                      double n_last) {
  int J = ctx->J;
  double total_n = P_prev * n_last;
  for (int j = 0; j < J; j++) {
    int idx_prev = j + (last - 1) * J;
    int idx_last = j + last * J;
    double base = ctx->cur_fpc_sum[idx_prev];
    double M_last = ctx->M[idx_last];
    double fpc = R_FINITE(M_last) ? (1.0 - n_last * ctx->rates[idx_last] / M_last) : 1.0;
    if (fpc < 0.0) fpc = 0.0;
    double cv2 = base + ctx->R[idx_last] * fpc / (total_n * ctx->rbar[idx_last]);
    if (cv_over(cv2, ctx->cv2_target[j])) return 0;
  }
  return 1;
}

static void eval_leaf(StageCtx *ctx, double P_prev, double C_prev) {
  int S = ctx->S;
  int J = ctx->J;
  int last = S - 1;
  double n_last;
  double lo = ctx->lower[last];
  double hi = ctx->upper[last];
  double c_last = ctx->cost[last];

  if (ctx->budget_mode) {
    double cap = budget_cap(ctx);
    double rem = cap - C_prev;
    if (rem < P_prev * c_last * lo) return;
    n_last = floor(rem / (P_prev * c_last));
    if (!R_FINITE(n_last)) return;
    if (R_FINITE(hi) && n_last > hi) n_last = hi;
    if (n_last < lo || n_last > INT_MAX_D) return;
    double total_n = P_prev * n_last;
    if (total_n < ctx->min_total - STAGE_TOL) return;
    double total_cost = C_prev + total_n * c_last;
    if (total_cost > cap) return;

    double max_ratio = 0.0;
    for (int j = 0; j < J; j++) {
      int idx_prev = j + (last - 1) * J;
      int idx_last = j + last * J;
      double base = ctx->cur_fpc_sum[idx_prev];
      double rbar_last = ctx->rbar[idx_last];
      double M_last = ctx->M[idx_last];
      double fpc = R_FINITE(M_last) ? (1.0 - n_last * ctx->rates[idx_last] / M_last) : 1.0;
      if (fpc < 0.0) fpc = 0.0;
      double cv2 = base + ctx->R[idx_last] * fpc / (total_n * rbar_last);
      double rat = cv2 / ctx->cv2_target[j];
      if (rat > max_ratio) max_ratio = rat;
    }

    double take_sum = n_last;
    for (int s = 1; s < last; s++) take_sum += ctx->cur_n[s];
    if (better(ctx, max_ratio, total_n, take_sum)) {
      store_incumbent(ctx, max_ratio, total_n, take_sum, total_cost,
                      max_ratio, n_last);
    }
  } else {
    /* CV mode: solve closed-form smallest feasible whole n_last */
    double need = lo;
    if (ctx->min_total > 0.0) {
      double mneed = ceil(ctx->min_total / P_prev - STAGE_TOL);
      if (mneed > need) need = mneed;
    }
    double floor_need = need;
    for (int j = 0; j < J; j++) {
      int idx_prev = j + (last - 1) * J;
      int idx_last = j + last * J;
      double base = ctx->cur_fpc_sum[idx_prev];
      double M_last = ctx->M[idx_last];
      double R_last = ctx->R[idx_last];
      double rbar_prev = ctx->rbar[idx_prev];
      double fpc_sub = R_FINITE(M_last) ? (R_last / (P_prev * rbar_prev * M_last)) : 0.0;
      double room = ctx->cv2_target[j] - (base - fpc_sub);
      if (R_last > 0.0) {
        if (room <= 0.0) return;
        double nj = R_last / (P_prev * ctx->rbar[idx_last] * room);
        if (!R_FINITE(nj)) return;
        double n_req = ceil(nj - STAGE_TOL);
        if (n_req > need) need = n_req;
      } else {
        if (cv_over(base, ctx->cv2_target[j])) return;
      }
    }
    n_last = need;
    if (!R_FINITE(n_last)) return;
    if ((R_FINITE(hi) && n_last > hi) || n_last > INT_MAX_D) return;
    if (!leaf_cv_ok(ctx, last, P_prev, n_last)) return;
    /* The closed form divides by a room one rounding error short of the
     * accepted one, which can amplify into a whole unit, so step down while
     * the exact check still holds. */
    while (n_last - 1.0 >= floor_need &&
           leaf_cv_ok(ctx, last, P_prev, n_last - 1.0)) {
      n_last -= 1.0;
    }
    double total_n = P_prev * n_last;
    double total_cost = C_prev + total_n * c_last;
    if (ctx->found && total_cost > cost_cap(ctx)) return;

    double take_sum = n_last;
    for (int s = 1; s < last; s++) take_sum += ctx->cur_n[s];
    if (better(ctx, total_cost, total_n, take_sum)) {
      store_incumbent(ctx, total_cost, total_n, take_sum, total_cost, 1.0,
                      n_last);
    }
  }
}

/* Evaluate a complete candidate prefix n[0..S-2] directly (for warm start / polish) */
static void eval_prefix(StageCtx *ctx, const double *prefix) {
  int S = ctx->S;
  int J = ctx->J;
  double P = 1.0;
  double C = 0.0;
  for (int s = 0; s < S - 1; s++) {
    double ns = prefix[s];
    if (!R_FINITE(ns)) return;
    if (ns < ctx->lower[s] || (R_FINITE(ctx->upper[s]) && ns > ctx->upper[s]) || ns > INT_MAX_D) {
      return;
    }
    ctx->cur_n[s] = ns;
    P *= ns;
    C += ctx->cost[s] * P;
    ctx->cur_P[s] = P;
    ctx->cur_C[s] = C;
    for (int j = 0; j < J; j++) {
      int idx = j + s * J;
      double prev_fpc = (s == 0) ? 0.0 : ctx->cur_fpc_sum[j + (s - 1) * J];
      double Ms = ctx->M[idx];
      double fpc = R_FINITE(Ms) ? (1.0 - ns * ctx->rates[idx] / Ms) : 1.0;
      if (fpc < 0.0) fpc = 0.0;
      ctx->cur_fpc_sum[idx] = prev_fpc + ctx->R[idx] * fpc / (P * ctx->rbar[idx]);
    }
  }
  eval_leaf(ctx, P, C);
}

/* Evaluate a single candidate value ns at stage s, updating state and recursing */
static void branch_value(StageCtx *ctx, int s, double ns, double P_prev, double C_prev) {
  int S = ctx->S;
  int J = ctx->J;
  double P_cur = P_prev * ns;
  double C_cur = C_prev + ctx->cost[s] * P_cur;

  /* Quick cost lower bound with remaining stages at their lower bounds */
  double min_complete_cost = C_prev + P_cur * ctx->rem_cost_min[s];
  if (ctx->budget_mode) {
    if (min_complete_cost > budget_cap(ctx)) return;
  } else {
    if (ctx->found && min_complete_cost > cost_cap(ctx)) return;
  }

  /* Update partial sums and check Cauchy-Schwarz & feasibility bounds */
  double rem_budget = budget_cap(ctx) - C_cur;
  for (int j = 0; j < J; j++) {
    int idx = j + s * J;
    double prev_fpc = (s == 0) ? 0.0 : ctx->cur_fpc_sum[j + (s - 1) * J];
    double prev_tel = (s == 0) ? (-ctx->kappa[j]) : ctx->cur_tel_sum[j + (s - 1) * J];
    double Ms = ctx->M[idx];
    double fpc = R_FINITE(Ms) ? (1.0 - ns * ctx->rates[idx] / Ms) : 1.0;
    if (fpc < 0.0) fpc = 0.0;
    double fpc_sum = prev_fpc + ctx->R[idx] * fpc / (P_cur * ctx->rbar[idx]);
    double tel_sum = prev_tel + ctx->alpha[idx] / P_cur;
    ctx->cur_fpc_sum[idx] = fpc_sum;
    ctx->cur_tel_sum[idx] = tel_sum;

    /* Check monotone lower bound at upper[s+1..S-1] */
    double min_cv2 = fpc_sum;
    double tail = ctx->g_tail[idx];
    if (tail > 0.0) {
      /* Subtract the stage-s contribution already counted in fpc_sum:
         g_tail[idx] holds sum_{k=s+1}^{S-1} R[j,k]*(1 - u_k*r/M)/(rbar[j,k]*prod_{m=s+1}^k u_m) */
      min_cv2 += tail / P_cur;
    }

    if (!ctx->budget_mode) {
      if (cv_over(min_cv2, ctx->cv2_target[j])) return;
      if (ctx->suffix_pos[idx]) {
        double room = ctx->cv2_target[j] * (1.0 + 1e-9) - tel_sum;
        if (room <= 0.0) {
          if (ctx->W[idx] > 0.0 || cv_over(tel_sum, ctx->cv2_target[j])) {
            return;
          }
        } else if (ctx->found && ctx->W[idx] > 0.0) {
          double cost_lb = C_cur + ctx->W[idx] / room;
          if (cost_lb > cost_cap(ctx)) return;
        }
      }
    } else {
      if (ctx->found) {
        double rat_lb = min_cv2 / ctx->cv2_target[j];
        if (rat_lb > ratio_cap(ctx)) return;
        if (ctx->suffix_pos[idx] && ctx->W[idx] > 0.0) {
          if (rem_budget <= 0.0) return;
          double cs_cv2 = tel_sum + ctx->W[idx] / rem_budget;
          double cs_rat = cs_cv2 / ctx->cv2_target[j];
          if (cs_rat > ratio_cap(ctx)) return;
        }
      }
    }
  }

  ctx->cur_n[s] = ns;
  ctx->cur_P[s] = P_cur;
  ctx->cur_C[s] = C_cur;

  if (s == S - 2) {
    eval_leaf(ctx, P_cur, C_cur);
    return;
  }

  /* Recurse to stage s + 1 */
  int next = s + 1;
  double lo = ctx->lower[next];
  double hi = ctx->upper[next];

  /* Cost upper bound on n_{next} */
  double cap = ctx->budget_mode ? budget_cap(ctx)
                                : (ctx->found ? cost_cap(ctx) : R_PosInf);
  if (R_FINITE(cap)) {
    double avail = cap - C_cur;
    if (avail < 0.0) return;
    double hi_cost = floor(avail / (P_cur * ctx->rem_cost_min[next]));
    if (!R_FINITE(hi) || hi_cost < hi) hi = hi_cost;
  }
  if (!R_FINITE(hi) || hi > INT_MAX_D) hi = INT_MAX_D;

  /* Variance lower bound on n_{next} in CV mode */
  if (!ctx->budget_mode) {
    for (int j = 0; j < J; j++) {
      int idx_cur = j + s * J;
      int idx_nxt = j + next * J;
      double prev_fpc = ctx->cur_fpc_sum[idx_cur];
      double Mn = ctx->M[idx_nxt];
      double Rn = ctx->R[idx_nxt];
      double sub = R_FINITE(Mn) ? (Rn / (P_cur * ctx->rbar[idx_cur] * Mn)) : 0.0;
      double g_num = Rn / ctx->rbar[idx_nxt] + ctx->g_tail[idx_nxt];
      /* The slack cv_over() accepts, so the bound never exceeds a size the
       * leaf would accept. */
      double room = ctx->cv2_target[j] * (1.0 + 1e-9) - (prev_fpc - sub);
      if (g_num > 0.0) {
        if (room <= 0.0) return;
        double req = ceil(g_num / (P_cur * room) - STAGE_TOL);
        if (req > lo) lo = req;
      } else {
        if (cv_over(prev_fpc, ctx->cv2_target[j])) return;
      }
    }
  }

  if (lo > hi) return;

  /* Scan outward from warm-start center for instant tight incumbent */
  double center = lo;
  if (ctx->n_init != NULL && R_FINITE(ctx->n_init[next])) {
    center = floor(ctx->n_init[next] + 0.5);
    if (center < lo) center = lo;
    if (center > hi) center = hi;
  }

  /* Visit center first, then center + d and center - d alternating */
  double max_d = (hi - center > center - lo) ? (hi - center) : (center - lo);
  for (double d = 0.0; d <= max_d; d += 1.0) {
    if (ctx->hit_limit) return;

    /* Dynamically tighten hi in CV mode if best_cost improved */
    if (!ctx->budget_mode && ctx->found) {
      double avail = cost_cap(ctx) - C_cur;
      if (avail < 0.0) return;
      double hi_cost = floor(avail / (P_cur * ctx->rem_cost_min[next]));
      if (hi_cost < hi) {
        hi = hi_cost;
        if (lo > hi) return;
        max_d = (hi - center > center - lo) ? (hi - center) : (center - lo);
        if (d > max_d) break;
      }
    }

    if (d == 0.0) {
      if (center >= lo && center <= hi) {
        ctx->nodes += 1.0;
        if (ctx->nodes > ctx->max_nodes) {
          ctx->hit_limit = 1;
          return;
        }
        if (((unsigned long long) ctx->nodes & 65535ULL) == 0ULL) R_CheckUserInterrupt();
        branch_value(ctx, next, center, P_cur, C_cur);
      }
    } else {
      double v_up = center + d;
      if (v_up <= hi && v_up >= lo) {
        ctx->nodes += 1.0;
        if (ctx->nodes > ctx->max_nodes) {
          ctx->hit_limit = 1;
          return;
        }
        if (((unsigned long long) ctx->nodes & 65535ULL) == 0ULL) R_CheckUserInterrupt();
        branch_value(ctx, next, v_up, P_cur, C_cur);
        if (ctx->hit_limit) return;
      }
      double v_dn = center - d;
      if (v_dn >= lo && v_dn <= hi) {
        ctx->nodes += 1.0;
        if (ctx->nodes > ctx->max_nodes) {
          ctx->hit_limit = 1;
          return;
        }
        if (((unsigned long long) ctx->nodes & 65535ULL) == 0ULL) R_CheckUserInterrupt();
        branch_value(ctx, next, v_dn, P_cur, C_cur);
      }
    }
  }
}

/* Warm start: try rounded corners of n_init and coordinate local polish */
static void warm_start(StageCtx *ctx) {
  int S = ctx->S;
  if (ctx->n_init == NULL) return;

  double *cand = (double *) R_alloc((size_t)(S - 1), sizeof(double));

  /* Try floor/ceil combinations of first S-1 stages (capped at S <= 10 -> 512 combos) */
  int n_free = S - 1;
  if (n_free <= 9) {
    int n_combos = 1 << n_free;
    for (int mask = 0; mask < n_combos; mask++) {
      for (int s = 0; s < n_free; s++) {
        double v = ctx->n_init[s];
        if (!R_FINITE(v)) v = ctx->lower[s];
        double vv = (mask & (1 << s)) ? ceil(v) : floor(v);
        if (vv < ctx->lower[s]) vv = ctx->lower[s];
        if (R_FINITE(ctx->upper[s]) && vv > ctx->upper[s]) vv = ctx->upper[s];
        cand[s] = vv;
      }
      eval_prefix(ctx, cand);
    }
  } else {
    for (int s = 0; s < n_free; s++) {
      double v = ctx->n_init[s];
      v = R_FINITE(v) ? floor(v + 0.5) : ctx->lower[s];
      if (v < ctx->lower[s]) v = ctx->lower[s];
      if (R_FINITE(ctx->upper[s]) && v > ctx->upper[s]) v = ctx->upper[s];
      cand[s] = v;
    }
    eval_prefix(ctx, cand);
  }

  /* Also try ceiling all takes upward (often feasible in CV mode) */
  for (int bump = 1; bump <= 4; bump++) {
    for (int s = 0; s < n_free; s++) {
      double v = ctx->n_init[s];
      v = (R_FINITE(v) ? ceil(v) : ctx->lower[s]) + (s == 0 ? bump : 0);
      if (v < ctx->lower[s]) v = ctx->lower[s];
      if (R_FINITE(ctx->upper[s]) && v > ctx->upper[s]) v = ctx->upper[s];
      cand[s] = v;
    }
    eval_prefix(ctx, cand);
  }

  /* Try upper bounds if finite (e.g. full census in budget or tight CV mode) */
  int all_fin = 1;
  for (int s = 0; s < n_free; s++) {
    if (!R_FINITE(ctx->upper[s])) {
      all_fin = 0;
      break;
    }
    cand[s] = ctx->upper[s];
  }
  if (all_fin) eval_prefix(ctx, cand);

  /* Coordinate local polish around best incumbent if found */
  if (ctx->found) {
    int improved = 1;
    int passes = 0;
    while (improved && passes < 50) {
      improved = 0;
      passes++;
      for (int s = 0; s < n_free; s++) {
        if (ctx->lower[s] >= ctx->upper[s]) continue;
        double base_v = ctx->best_n[s];
        int steps[8] = {-1, 1, -2, 2, -4, 4, -8, 8};
        for (int k = 0; k < 8; k++) {
          double v = base_v + steps[k];
          if (v < ctx->lower[s]) continue;
          if (R_FINITE(ctx->upper[s]) && v > ctx->upper[s]) continue;
          for (int q = 0; q < n_free; q++) cand[q] = ctx->best_n[q];
          cand[s] = v;
          long before = ctx->n_improve;
          eval_prefix(ctx, cand);
          if (ctx->n_improve != before) {
            improved = 1;
            base_v = ctx->best_n[s];
          }
        }
      }
    }
  }
}

static void check_vector(SEXP x, const char *name, int len) {
  if (!Rf_isReal(x) || XLENGTH(x) != len) {
    Rf_error("'%s' must be a double vector of length %d", name, len);
  }
}

static void check_matrix(SEXP x, const char *name, int J, int S) {
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  if (!Rf_isReal(x) || Rf_isNull(d) || Rf_length(d) != 2 ||
      INTEGER(d)[0] != J || INTEGER(d)[1] != S) {
    Rf_error("'%s' must be a %d x %d double matrix", name, J, S);
  }
}

SEXP C_stage_search(
  SEXP R_sexp,
  SEXP rates_sexp,
  SEXP M_sexp,
  SEXP cv2_target_sexp,
  SEXP cost_sexp,
  SEXP lower_sexp,
  SEXP upper_sexp,
  SEXP budget_var_sexp,
  SEXP min_total_sexp,
  SEXP n_init_sexp,
  SEXP max_nodes_sexp
) {
  SEXP dims = Rf_getAttrib(R_sexp, R_DimSymbol);
  if (!Rf_isReal(R_sexp) || Rf_isNull(dims) || Rf_length(dims) != 2) {
    Rf_error("R must be a J x S double matrix");
  }
  int J = INTEGER(dims)[0];
  int S = INTEGER(dims)[1];
  if (J < 1 || S < 2) {
    Rf_error("R must have at least 1 row and 2 columns");
  }
  check_matrix(rates_sexp, "rates", J, S);
  check_matrix(M_sexp, "M", J, S);
  check_vector(cv2_target_sexp, "cv2_target", J);
  check_vector(cost_sexp, "cost", S);
  check_vector(lower_sexp, "lower", S);
  check_vector(upper_sexp, "upper", S);
  if (!Rf_isNull(n_init_sexp)) check_vector(n_init_sexp, "n_init", S);
  for (int s = 0; s < S; s++) {
    double lo = REAL(lower_sexp)[s];
    double hi = REAL(upper_sexp)[s];
    if (!R_FINITE(lo) || lo < 1.0 || ISNAN(hi)) {
      Rf_error("'lower' must be finite and at least 1, and 'upper' not NaN");
    }
  }

  StageCtx ctx;
  ctx.S = S;
  ctx.J = J;
  ctx.budget_var = Rf_asReal(budget_var_sexp);
  ctx.budget_mode = !ISNA(ctx.budget_var) && R_FINITE(ctx.budget_var);
  ctx.min_total = Rf_asReal(min_total_sexp);
  ctx.max_nodes = Rf_asReal(max_nodes_sexp);
  if (!R_FINITE(ctx.max_nodes) || ctx.max_nodes <= 0.0) ctx.max_nodes = 2e7;
  ctx.nodes = 0.0;
  ctx.hit_limit = 0;

  ctx.R = REAL(R_sexp);
  ctx.rates = REAL(rates_sexp);
  ctx.M = REAL(M_sexp);
  ctx.cv2_target = REAL(cv2_target_sexp);
  ctx.cost = REAL(cost_sexp);
  ctx.lower = REAL(lower_sexp);
  ctx.upper = REAL(upper_sexp);
  ctx.n_init = Rf_isNull(n_init_sexp) ? NULL : REAL(n_init_sexp);

  ctx.rbar = (double *) R_alloc((size_t)(J * S), sizeof(double));
  ctx.alpha = (double *) R_alloc((size_t)(J * S), sizeof(double));
  ctx.kappa = (double *) R_alloc((size_t)J, sizeof(double));
  ctx.W = (double *) R_alloc((size_t)(J * S), sizeof(double));
  ctx.suffix_pos = (int *) R_alloc((size_t)(J * S), sizeof(int));
  ctx.g_tail = (double *) R_alloc((size_t)(J * S), sizeof(double));
  ctx.rem_cost_min = (double *) R_alloc((size_t)S, sizeof(double));

  ctx.cur_n = (double *) R_alloc((size_t)S, sizeof(double));
  ctx.cur_P = (double *) R_alloc((size_t)S, sizeof(double));
  ctx.cur_C = (double *) R_alloc((size_t)S, sizeof(double));
  ctx.cur_fpc_sum = (double *) R_alloc((size_t)(J * S), sizeof(double));
  ctx.cur_tel_sum = (double *) R_alloc((size_t)(J * S), sizeof(double));
  ctx.best_n = (double *) R_alloc((size_t)S, sizeof(double));

  ctx.found = 0;
  ctx.n_improve = 0;
  ctx.best_prim = R_PosInf;
  ctx.best_tn = R_PosInf;
  ctx.best_take = R_PosInf;
  ctx.best_cost = R_PosInf;
  ctx.best_ratio = R_PosInf;

  /* Precompute cumulative rates, telescoping alpha, kappa */
  for (int j = 0; j < J; j++) {
    double rb = 1.0;
    for (int s = 0; s < S; s++) {
      int idx = j + s * J;
      rb *= ctx.rates[idx];
      ctx.rbar[idx] = rb;
    }
    ctx.kappa[j] = R_FINITE(ctx.M[j]) ? (ctx.R[j] / ctx.M[j]) : 0.0;
    for (int s = 0; s < S; s++) {
      int idx = j + s * J;
      double r_cur = ctx.R[idx];
      double sub = 0.0;
      if (s < S - 1) {
        int idx_nxt = j + (s + 1) * J;
        if (R_FINITE(ctx.M[idx_nxt])) {
          sub = ctx.R[idx_nxt] / ctx.M[idx_nxt];
        }
      }
      ctx.alpha[idx] = (r_cur - sub) / ctx.rbar[idx];
    }
  }

  /* Precompute Cauchy-Schwarz suffix weights W[j, s] for k = s+1 .. S-1 */
  for (int j = 0; j < J; j++) {
    for (int s = 0; s < S; s++) {
      int idx = j + s * J;
      int all_pos = 1;
      double sum_sq = 0.0;
      for (int k = s + 1; k < S; k++) {
        double a = ctx.alpha[j + k * J];
        if (a < 0.0) {
          all_pos = 0;
          break;
        }
        sum_sq += sqrt(a * ctx.cost[k]);
      }
      ctx.suffix_pos[idx] = all_pos;
      ctx.W[idx] = all_pos ? (sum_sq * sum_sq) : 0.0;
    }
  }

  /* Precompute monotone variance tails g_tail[j, s] evaluated at upper[s+1..S-1] */
  for (int j = 0; j < J; j++) {
    ctx.g_tail[j + (S - 1) * J] = 0.0;
    for (int s = S - 2; s >= 0; s--) {
      double tail = 0.0;
      double prod_u = 1.0;
      for (int k = s + 1; k < S; k++) {
        double uk = ctx.upper[k];
        if (!R_FINITE(uk)) break;
        prod_u *= uk;
        int idx_k = j + k * J;
        double Mk = ctx.M[idx_k];
        double fpc = R_FINITE(Mk) ? (1.0 - uk * ctx.rates[idx_k] / Mk) : 1.0;
        if (fpc < 0.0) fpc = 0.0;
        tail += ctx.R[idx_k] * fpc / (ctx.rbar[idx_k] * prod_u);
      }
      ctx.g_tail[j + s * J] = tail;
    }
  }

  /* Precompute minimum remaining cost per unit of P_s:
     rem_cost_min[s] = cost[s] + lower[s+1] * rem_cost_min[s+1] */
  ctx.rem_cost_min[S - 1] = ctx.cost[S - 1];
  for (int s = S - 2; s >= 0; s--) {
    ctx.rem_cost_min[s] = ctx.cost[s] + ctx.lower[s + 1] * ctx.rem_cost_min[s + 1];
  }

  /* Warm start from continuous optimum */
  warm_start(&ctx);

  /* Stage-0 bounds */
  double lo0 = ctx.lower[0];
  double hi0 = ctx.upper[0];
  if (ctx.budget_mode) {
    double hi_cost = floor(budget_cap(&ctx) / ctx.rem_cost_min[0]);
    if (!R_FINITE(hi0) || hi_cost < hi0) hi0 = hi_cost;
  } else {
    for (int j = 0; j < J; j++) {
      double g0 = ctx.R[j] / ctx.rbar[j] + ctx.g_tail[j];
      double denom = ctx.cv2_target[j] * (1.0 + 1e-9) + ctx.kappa[j];
      if (g0 > 0.0) {
        if (denom <= 0.0) {
          lo0 = R_PosInf;
          break;
        }
        double req = ceil(g0 / denom - STAGE_TOL);
        if (req > lo0) lo0 = req;
      }
    }
    if (ctx.found) {
      double hi_cost = floor(cost_cap(&ctx) / ctx.rem_cost_min[0]);
      if (!R_FINITE(hi0) || hi_cost < hi0) hi0 = hi_cost;
    }
  }
  if (!R_FINITE(hi0) || hi0 > INT_MAX_D) hi0 = INT_MAX_D;

  /* In CV mode without warm start, seed a feasible stage-0 count first so hi0 becomes finite */
  if (!ctx.budget_mode && !ctx.found && R_FINITE(lo0) && lo0 <= hi0) {
    double *seed = (double *) R_alloc((size_t)(S - 1), sizeof(double));
    double a_seed = lo0;
    for (int s = 1; s < S - 1; s++) {
      seed[s] = (ctx.n_init != NULL && R_FINITE(ctx.n_init[s]))
        ? floor(ctx.n_init[s] + 0.5)
        : ctx.lower[s];
      if (seed[s] < ctx.lower[s]) seed[s] = ctx.lower[s];
      if (R_FINITE(ctx.upper[s]) && seed[s] > ctx.upper[s]) seed[s] = ctx.upper[s];
    }
    /* Find a stage-0 count above alpha[0] + ... at seed takes */
    for (int j = 0; j < J; j++) {
      double g_take = ctx.alpha[j];
      double pt = 1.0;
      for (int s = 1; s < S - 1; s++) {
        pt *= seed[s];
        g_take += ctx.alpha[j + s * J] / pt;
      }
      double req = floor(g_take / (ctx.cv2_target[j] + ctx.kappa[j]) + STAGE_TOL) + 1.0;
      if (req > a_seed) a_seed = req;
    }
    if (R_FINITE(ctx.upper[0]) && a_seed > ctx.upper[0]) a_seed = ctx.upper[0];
    if (a_seed >= lo0 && a_seed <= hi0) {
      seed[0] = a_seed;
      eval_prefix(&ctx, seed);
      if (ctx.found) {
        double hi_cost = floor(cost_cap(&ctx) / ctx.rem_cost_min[0]);
        if (hi_cost < hi0) hi0 = hi_cost;
      }
    }
  }

  /* Branch-and-bound over stage 0, scanning outward from n_init[0] */
  if (R_FINITE(lo0) && lo0 <= hi0) {
    double center = lo0;
    if (ctx.n_init != NULL && R_FINITE(ctx.n_init[0])) {
      center = floor(ctx.n_init[0] + 0.5);
      if (center < lo0) center = lo0;
      if (center > hi0) center = hi0;
    }
    double max_d = (hi0 - center > center - lo0) ? (hi0 - center) : (center - lo0);
    for (double d = 0.0; d <= max_d; d += 1.0) {
      if (ctx.hit_limit) break;

      if (!ctx.budget_mode && ctx.found) {
        double hi_cost = floor(cost_cap(&ctx) / ctx.rem_cost_min[0]);
        if (hi_cost < hi0) {
          hi0 = hi_cost;
          if (lo0 > hi0) break;
          max_d = (hi0 - center > center - lo0) ? (hi0 - center) : (center - lo0);
          if (d > max_d) break;
        }
      }

      if (d == 0.0) {
        if (center >= lo0 && center <= hi0) {
          ctx.nodes += 1.0;
          if (ctx.nodes > ctx.max_nodes) {
            ctx.hit_limit = 1;
            break;
          }
          branch_value(&ctx, 0, center, 1.0, 0.0);
        }
      } else {
        double v_up = center + d;
        if (v_up <= hi0 && v_up >= lo0) {
          ctx.nodes += 1.0;
          if (ctx.nodes > ctx.max_nodes) {
            ctx.hit_limit = 1;
            break;
          }
          if (((unsigned long long) ctx.nodes & 65535ULL) == 0ULL) R_CheckUserInterrupt();
          branch_value(&ctx, 0, v_up, 1.0, 0.0);
          if (ctx.hit_limit) break;
        }
        double v_dn = center - d;
        if (v_dn >= lo0 && v_dn <= hi0) {
          ctx.nodes += 1.0;
          if (ctx.nodes > ctx.max_nodes) {
            ctx.hit_limit = 1;
            break;
          }
          if (((unsigned long long) ctx.nodes & 65535ULL) == 0ULL) R_CheckUserInterrupt();
          branch_value(&ctx, 0, v_dn, 1.0, 0.0);
        }
      }
    }
  }

  if (ctx.hit_limit && ctx.nodes > ctx.max_nodes) ctx.nodes = ctx.max_nodes;

  /* Package result list: n (integer S), found (logical), exact (logical), nodes (double) */
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 4));
  SET_STRING_ELT(nms, 0, Rf_mkChar("n"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("found"));
  SET_STRING_ELT(nms, 2, Rf_mkChar("exact"));
  SET_STRING_ELT(nms, 3, Rf_mkChar("nodes"));
  Rf_setAttrib(out, R_NamesSymbol, nms);

  SEXP n_out = PROTECT(Rf_allocVector(INTSXP, S));
  int *np = INTEGER(n_out);
  if (ctx.found) {
    for (int s = 0; s < S; s++) np[s] = (int) ctx.best_n[s];
  } else {
    for (int s = 0; s < S; s++) np[s] = NA_INTEGER;
  }
  SET_VECTOR_ELT(out, 0, n_out);
  SET_VECTOR_ELT(out, 1, Rf_ScalarLogical(ctx.found));
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(!ctx.hit_limit));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(ctx.nodes));

  UNPROTECT(3);
  return out;
}
