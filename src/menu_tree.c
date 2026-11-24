/* Exact branch-and-bound solver for the group-selection menu problem.
 *
 * Implements .menu_tree() in compiled C with:
 * - Row-major constraint storage for cache-contiguous dot products.
 * - Pre-sorted per-group states by cost for O(1)-break leaf scans.
 * - Binary min-heap for best-first open nodes and LIFO stack for dives.
 * - Fused single-pass Lagrangian subgradient and exact 1D piecewise-linear
 *   dual maximization without intermediate heap allocations.
 * - Interrupt-safe memory management via R_ToplevelExec().
 */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  double L;
  double c_fix;
  int branch;
  int total_alive;
  int *plan;       /* length H: -1 if free, else compacted state index */
  double *cap;     /* length K */
  double *mu;      /* length K */
  int *alive_off;  /* length H + 1 */
  int *alive_idx;  /* length total_alive */
} MenuNode;

typedef struct {
  int H;
  int K;
  int N;           /* compacted kept states */
  int kept_total;
  double *cost;    /* length N */
  double *C_row;   /* N * K row-major */
  int *orig_idx;   /* length N: 1-based original state index */
  int *grp_off;    /* length H + 1 */

  double grain;
  double tol;
  double ceiling_cost;
  double best;
  int *best_plan;  /* length H: compacted state index, or -1 */
  int has_plan;
  int exact;

  /* Open min-heap and DFS stack */
  MenuNode **open;
  int open_size;
  int open_cap;
  MenuNode **stack;
  int stack_size;
  int stack_cap;

  /* Scratch buffers for evaluate_node (reused across all nodes) */
  int *free_list;  /* length H */
  double *cmin;    /* length K */
  double *mu_cur;  /* length K */
  double *mu_best; /* length K */
  double *grad;    /* length K */
  int *grp_arg;    /* length H */
  int *tmp_off;    /* length H + 1 */
  int *tmp_idx;    /* length N */
  int *tmp_plan;   /* length H */

  /* Children scratch */
  MenuNode **kids;
  int kids_cap;
} MenuCtx;

static void check_interrupt_cb(void *dummy) {
  (void) dummy;
  R_CheckUserInterrupt();
}

static int check_interrupted(void) {
  return !R_ToplevelExec(check_interrupt_cb, NULL);
}

static void free_node(MenuNode *nd) {
  if (!nd) return;
  free(nd->plan);
  free(nd->cap);
  free(nd->mu);
  free(nd->alive_off);
  free(nd->alive_idx);
  free(nd);
}

static MenuNode *alloc_node(int H, int K, int total_alive) {
  MenuNode *nd = (MenuNode *) malloc(sizeof(MenuNode));
  if (!nd) return NULL;
  nd->plan = (int *) malloc((size_t) H * sizeof(int));
  nd->cap = (double *) malloc((size_t) (K > 0 ? K : 1) * sizeof(double));
  nd->mu = (double *) malloc((size_t) (K > 0 ? K : 1) * sizeof(double));
  nd->alive_off = (int *) malloc((size_t) (H + 1) * sizeof(int));
  nd->alive_idx = (int *) malloc((size_t) (total_alive > 0 ? total_alive : 1) * sizeof(int));
  if (!nd->plan || !nd->cap || !nd->mu || !nd->alive_off || !nd->alive_idx) {
    free_node(nd);
    return NULL;
  }
  nd->L = R_NegInf;
  nd->c_fix = 0.0;
  nd->branch = -1;
  nd->total_alive = total_alive;
  return nd;
}

static void free_ctx(MenuCtx *ctx) {
  for (int i = 0; i < ctx->open_size; i++) free_node(ctx->open[i]);
  for (int i = 0; i < ctx->stack_size; i++) free_node(ctx->stack[i]);
  free(ctx->open);
  free(ctx->stack);
  free(ctx->kids);
  free(ctx->cost);
  free(ctx->C_row);
  free(ctx->orig_idx);
  free(ctx->grp_off);
  free(ctx->best_plan);
  free(ctx->free_list);
  free(ctx->cmin);
  free(ctx->mu_cur);
  free(ctx->mu_best);
  free(ctx->grad);
  free(ctx->grp_arg);
  free(ctx->tmp_off);
  free(ctx->tmp_idx);
  free(ctx->tmp_plan);
}

static inline double ctx_limit(const MenuCtx *ctx) {
  if (R_FINITE(ctx->best)) {
    double step = ctx->grain - ctx->tol;
    if (step < ctx->tol) step = ctx->tol;
    return ctx->best - step;
  }
  return ctx->ceiling_cost;
}

static inline void ctx_improve(MenuCtx *ctx, double val, const int *plan) {
  if (val < ctx->best - ctx->tol) {
    ctx->best = val;
    memcpy(ctx->best_plan, plan, (size_t) ctx->H * sizeof(int));
    ctx->has_plan = 1;
    double scale = fabs(val);
    if (scale < 1e-300) scale = 1e-300;
    ctx->tol = 1e-9 * scale;
  }
}

/* Min-heap helpers on ctx->open keyed by nd->L */
static int heap_push(MenuCtx *ctx, MenuNode *nd) {
  if (ctx->open_size >= ctx->open_cap) {
    int ncap = ctx->open_cap == 0 ? 64 : ctx->open_cap * 2;
    MenuNode **tmp = (MenuNode **) realloc(ctx->open, (size_t) ncap * sizeof(MenuNode *));
    if (!tmp) return 0;
    ctx->open = tmp;
    ctx->open_cap = ncap;
  }
  int i = ctx->open_size++;
  while (i > 0) {
    int p = (i - 1) >> 1;
    if (ctx->open[p]->L <= nd->L) break;
    ctx->open[i] = ctx->open[p];
    i = p;
  }
  ctx->open[i] = nd;
  return 1;
}

static MenuNode *heap_pop(MenuCtx *ctx) {
  MenuNode *top = ctx->open[0];
  int n = --ctx->open_size;
  if (n > 0) {
    MenuNode *last = ctx->open[n];
    int i = 0;
    for (;;) {
      int l = (i << 1) + 1;
      if (l >= n) break;
      int r = l + 1;
      int s = (r < n && ctx->open[r]->L < ctx->open[l]->L) ? r : l;
      if (last->L <= ctx->open[s]->L) break;
      ctx->open[i] = ctx->open[s];
      i = s;
    }
    ctx->open[i] = last;
  }
  return top;
}

static int stack_push(MenuCtx *ctx, MenuNode *nd) {
  if (ctx->stack_size >= ctx->stack_cap) {
    int ncap = ctx->stack_cap == 0 ? 64 : ctx->stack_cap * 2;
    MenuNode **tmp = (MenuNode **) realloc(ctx->stack, (size_t) ncap * sizeof(MenuNode *));
    if (!tmp) return 0;
    ctx->stack = tmp;
    ctx->stack_cap = ncap;
  }
  ctx->stack[ctx->stack_size++] = nd;
  return 1;
}

/* Sorting helpers for root dominance & cost order */
typedef struct {
  int orig;
  double cost;
  double csum;
} SortItem;

static int cmp_sort_item(const void *a, const void *b) {
  const SortItem *sa = (const SortItem *) a;
  const SortItem *sb = (const SortItem *) b;
  if (sa->cost < sb->cost) return -1;
  if (sa->cost > sb->cost) return 1;
  if (sa->csum < sb->csum) return -1;
  if (sa->csum > sb->csum) return 1;
  return (sa->orig < sb->orig) ? -1 : (sa->orig > sb->orig);
}

/* Children by descending bound, ties kept in creation order, so the dive
 * order and any result an early stop leaves do not depend on the C
 * library's qsort. */
static void sort_kids_desc(MenuNode **kids, int n) {
  for (int a = 1; a < n; a++) {
    MenuNode *key = kids[a];
    int b = a - 1;
    while (b >= 0 && kids[b]->L < key->L) {
      kids[b + 1] = kids[b];
      b--;
    }
    kids[b + 1] = key;
  }
}

/* Evaluate 1D Lagrangian dual at multiplier m for K == 1 */
static void dual_one_eval(MenuCtx *ctx, const MenuNode *nd, int n_free,
                          double m, double *out_L, double *out_g) {
  const double *cost = ctx->cost;
  const double *C_row = ctx->C_row;
  double sum_min = 0.0;
  double sum_c = 0.0;
  for (int fi = 0; fi < n_free; fi++) {
    int h = ctx->free_list[fi];
    int a0 = nd->alive_off[h];
    int a1 = nd->alive_off[h + 1];
    int s0 = nd->alive_idx[a0];
    double m1 = cost[s0] + C_row[s0] * m;
    int arg = s0;
    for (int p = a0 + 1; p < a1; p++) {
      int s = nd->alive_idx[p];
      double v = cost[s] + C_row[s] * m;
      if (v < m1) {
        m1 = v;
        arg = s;
      }
    }
    sum_min += m1;
    sum_c += C_row[arg];
    ctx->grp_arg[fi] = arg;
  }
  double g = sum_c - nd->cap[0];
  if (g <= 1e-9) {
    double total_c = nd->c_fix;
    memcpy(ctx->tmp_plan, nd->plan, (size_t) ctx->H * sizeof(int));
    for (int fi = 0; fi < n_free; fi++) {
      int arg = ctx->grp_arg[fi];
      ctx->tmp_plan[ctx->free_list[fi]] = arg;
      total_c += cost[arg];
    }
    ctx_improve(ctx, total_c, ctx->tmp_plan);
  }
  *out_L = nd->c_fix + sum_min - m * nd->cap[0];
  *out_g = g;
}

static int solve_dual_one(MenuCtx *ctx, const MenuNode *nd, int n_free,
                          double lambda0, double *best_L_out,
                          double *best_lam_out) {
  double fa_L, fa_g;
  dual_one_eval(ctx, nd, n_free, 0.0, &fa_L, &fa_g);
  if (fa_g <= 1e-9) {
    *best_L_out = fa_L;
    *best_lam_out = 0.0;
    return 1;
  }
  double a = 0.0;
  double b = lambda0 > 1e-300 ? lambda0 : 1e-300;
  double fb_L, fb_g;
  dual_one_eval(ctx, nd, n_free, b, &fb_L, &fb_g);
  int steps = 2;
  while (fb_g > 1e-9) {
    if (steps > 400) return 0;
    a = b;
    fa_L = fb_L;
    fa_g = fb_g;
    b *= 4.0;
    if (!R_FINITE(b)) return 0;
    dual_one_eval(ctx, nd, n_free, b, &fb_L, &fb_g);
    steps++;
  }
  double best_L, best_lam;
  if (fa_L >= fb_L) {
    best_L = fa_L;
    best_lam = a;
  } else {
    best_L = fb_L;
    best_lam = b;
  }
  for (int it = 0; it < 200; it++) {
    double denom = fa_g - fb_g;
    if (denom == 0.0) break;
    double m = (fb_L - fa_L + fa_g * a - fb_g * b) / denom;
    if (!(m > a && m < b)) break;
    double v = fa_L + fa_g * (m - a);
    double fm_L, fm_g;
    dual_one_eval(ctx, nd, n_free, m, &fm_L, &fm_g);
    if (fm_L > best_L) {
      best_L = fm_L;
      best_lam = m;
    }
    if (best_L > ctx_limit(ctx)) break;
    double av = fabs(v);
    if (v - fm_L <= 1e-12 * (av > 1e-300 ? av : 1e-300)) break;
    if (fm_g > 1e-9) {
      a = m;
      fa_L = fm_L;
      fa_g = fm_g;
    } else {
      b = m;
      fb_L = fm_L;
      fb_g = fm_g;
    }
  }
  *best_L_out = best_L;
  *best_lam_out = best_lam;
  return 1;
}

/* Evaluate a node in-place. Returns nd if it survives, or frees nd and returns NULL. */
static MenuNode *evaluate_node(MenuCtx *ctx, MenuNode *nd) {
  const int H = ctx->H;
  const int K = ctx->K;
  const double *cost = ctx->cost;
  const double *C_row = ctx->C_row;

  for (;;) {
    /* 1. Fix any free group reduced to 1 alive state */
    int fixed_any;
    do {
      fixed_any = 0;
      for (int h = 0; h < H; h++) {
        if (nd->plan[h] < 0 && (nd->alive_off[h + 1] - nd->alive_off[h]) == 1) {
          int s = nd->alive_idx[nd->alive_off[h]];
          nd->c_fix += cost[s];
          const double *crow = C_row + (size_t) s * K;
          for (int k = 0; k < K; k++) nd->cap[k] -= crow[k];
          nd->plan[h] = s;
          fixed_any = 1;
        }
      }
    } while (fixed_any);

    /* 2. Check residual capacity */
    for (int k = 0; k < K; k++) {
      if (nd->cap[k] < -1e-9) {
        free_node(nd);
        return NULL;
      }
    }

    /* 3. Collect free groups */
    int n_free = 0;
    for (int h = 0; h < H; h++) {
      if (nd->plan[h] < 0) ctx->free_list[n_free++] = h;
    }

    if (n_free == 0) {
      ctx_improve(ctx, nd->c_fix, nd->plan);
      free_node(nd);
      return NULL;
    }

    /* 4. Exact 1-free-group leaf (states already sorted by cost ascending!) */
    if (n_free == 1) {
      int h0 = ctx->free_list[0];
      int a0 = nd->alive_off[h0];
      int a1 = nd->alive_off[h0 + 1];
      for (int p = a0; p < a1; p++) {
        int s = nd->alive_idx[p];
        const double *crow = C_row + (size_t) s * K;
        int ok = 1;
        for (int k = 0; k < K; k++) {
          if (crow[k] > nd->cap[k] + 1e-9) {
            ok = 0;
            break;
          }
        }
        if (ok) {
          nd->plan[h0] = s;
          ctx_improve(ctx, nd->c_fix + cost[s], nd->plan);
          break;
        }
      }
      free_node(nd);
      return NULL;
    }

    /* 5. Exact 2-free-group leaf (both groups sorted by cost ascending!) */
    if (n_free == 2) {
      int g1 = ctx->free_list[0];
      int g2 = ctx->free_list[1];
      if ((nd->alive_off[g2 + 1] - nd->alive_off[g2]) <
          (nd->alive_off[g1 + 1] - nd->alive_off[g1])) {
        int tmp = g1; g1 = g2; g2 = tmp;
      }
      int a1_0 = nd->alive_off[g1], a1_1 = nd->alive_off[g1 + 1];
      int a2_0 = nd->alive_off[g2], a2_1 = nd->alive_off[g2 + 1];
      double min_c2 = cost[nd->alive_idx[a2_0]];
      for (int p1 = a1_0; p1 < a1_1; p1++) {
        int s1 = nd->alive_idx[p1];
        double c1 = nd->c_fix + cost[s1];
        if (c1 + min_c2 >= ctx->best - ctx->tol) break;
        const double *r1 = C_row + (size_t) s1 * K;
        for (int p2 = a2_0; p2 < a2_1; p2++) {
          int s2 = nd->alive_idx[p2];
          if (c1 + cost[s2] >= ctx->best - ctx->tol) break;
          const double *r2 = C_row + (size_t) s2 * K;
          int ok = 1;
          for (int k = 0; k < K; k++) {
            if (r1[k] + r2[k] > nd->cap[k] + 1e-9) {
              ok = 0;
              break;
            }
          }
          if (ok) {
            nd->plan[g1] = s1;
            nd->plan[g2] = s2;
            ctx_improve(ctx, c1 + cost[s2], nd->plan);
            break;
          }
        }
      }
      free_node(nd);
      return NULL;
    }

    /* 6. Lower-bound feasibility check via per-group column minima */
    for (int k = 0; k < K; k++) ctx->cmin[k] = 0.0;
    for (int fi = 0; fi < n_free; fi++) {
      int h = ctx->free_list[fi];
      int a0 = nd->alive_off[h];
      int a1 = nd->alive_off[h + 1];
      const double *r0 = C_row + (size_t) nd->alive_idx[a0] * K;
      for (int k = 0; k < K; k++) ctx->grad[k] = r0[k];
      for (int p = a0 + 1; p < a1; p++) {
        const double *rp = C_row + (size_t) nd->alive_idx[p] * K;
        for (int k = 0; k < K; k++) {
          if (rp[k] < ctx->grad[k]) ctx->grad[k] = rp[k];
        }
      }
      for (int k = 0; k < K; k++) ctx->cmin[k] += ctx->grad[k];
    }
    for (int k = 0; k < K; k++) {
      if (ctx->cmin[k] > nd->cap[k] + 1e-9) {
        free_node(nd);
        return NULL;
      }
    }

    /* 7. Lagrangian dual bound (exact 1D or 30 subgradient iterations) */
    double best_L = R_NegInf;
    int solved_exact = 0;
    if (ctx->exact && K == 1) {
      double c_lo = R_PosInf, c_hi = R_NegInf, max_cv = 1e-300;
      for (int fi = 0; fi < n_free; fi++) {
        int h = ctx->free_list[fi];
        int a0 = nd->alive_off[h];
        int a1 = nd->alive_off[h + 1];
        double c_first = cost[nd->alive_idx[a0]];
        double c_last = cost[nd->alive_idx[a1 - 1]];
        if (c_first < c_lo) c_lo = c_first;
        if (c_last > c_hi) c_hi = c_last;
        for (int p = a0; p < a1; p++) {
          double cv = C_row[nd->alive_idx[p]];
          if (cv > max_cv) max_cv = cv;
        }
      }
      double scale = (c_hi - c_lo) / max_cv;
      double lam0 = nd->mu[0];
      if (scale > lam0) lam0 = scale;
      if (1e-300 > lam0) lam0 = 1e-300;
      double lam_opt = 0.0;
      if (solve_dual_one(ctx, nd, n_free, lam0, &best_L, &lam_opt)) {
        ctx->mu_best[0] = lam_opt;
        solved_exact = 1;
      }
    }

    if (!solved_exact) {
      memcpy(ctx->mu_cur, nd->mu, (size_t) K * sizeof(double));
      memcpy(ctx->mu_best, nd->mu, (size_t) K * sizeof(double));
      double theta = 1.0;
      for (int it = 0; it < 30; it++) {
        double sum_min = 0.0;
        for (int k = 0; k < K; k++) ctx->grad[k] = -nd->cap[k];
        for (int fi = 0; fi < n_free; fi++) {
          int h = ctx->free_list[fi];
          int a0 = nd->alive_off[h];
          int a1 = nd->alive_off[h + 1];
          int s0 = nd->alive_idx[a0];
          const double *r0 = C_row + (size_t) s0 * K;
          double m1 = cost[s0];
          for (int k = 0; k < K; k++) m1 += r0[k] * ctx->mu_cur[k];
          int arg = s0;
          for (int p = a0 + 1; p < a1; p++) {
            int s = nd->alive_idx[p];
            const double *rp = C_row + (size_t) s * K;
            double v = cost[s];
            for (int k = 0; k < K; k++) v += rp[k] * ctx->mu_cur[k];
            if (v < m1) {
              m1 = v;
              arg = s;
            }
          }
          sum_min += m1;
          ctx->grp_arg[fi] = arg;
          const double *ra = C_row + (size_t) arg * K;
          for (int k = 0; k < K; k++) ctx->grad[k] += ra[k];
        }
        double mu_cap = 0.0;
        for (int k = 0; k < K; k++) mu_cap += ctx->mu_cur[k] * nd->cap[k];
        double L = nd->c_fix + sum_min - mu_cap;
        if (L > best_L) {
          best_L = L;
          memcpy(ctx->mu_best, ctx->mu_cur, (size_t) K * sizeof(double));
        } else {
          theta *= 0.7;
        }
        int feas = 1;
        for (int k = 0; k < K; k++) {
          if (ctx->grad[k] > 1e-9) {
            feas = 0;
            break;
          }
        }
        if (feas) {
          double total_c = nd->c_fix;
          memcpy(ctx->tmp_plan, nd->plan, (size_t) H * sizeof(int));
          for (int fi = 0; fi < n_free; fi++) {
            int arg = ctx->grp_arg[fi];
            ctx->tmp_plan[ctx->free_list[fi]] = arg;
            total_c += cost[arg];
          }
          ctx_improve(ctx, total_c, ctx->tmp_plan);
        }
        double lim = ctx_limit(ctx);
        if (best_L > lim) break;
        double g2 = 0.0;
        for (int k = 0; k < K; k++) {
          if (ctx->mu_cur[k] <= 0.0 && ctx->grad[k] < 0.0) ctx->grad[k] = 0.0;
          g2 += ctx->grad[k] * ctx->grad[k];
        }
        if (g2 == 0.0) break;
        double step = theta * (lim - L) / g2;
        for (int k = 0; k < K; k++) {
          double nm = ctx->mu_cur[k] + step * ctx->grad[k];
          ctx->mu_cur[k] = nm > 0.0 ? nm : 0.0;
        }
      }
    }

    double slack = ctx_limit(ctx) - best_L;
    if (slack < 0.0) {
      free_node(nd);
      return NULL;
    }

    /* 8. Reduced-cost filtering at best_mu + branch selection */
    int new_total = 0;
    int branch_h = -1;
    double best_gap = R_PosInf;
    for (int h = 0; h < H; h++) {
      ctx->tmp_off[h] = new_total;
      if (nd->plan[h] >= 0) continue;
      int a0 = nd->alive_off[h];
      int a1 = nd->alive_off[h + 1];
      double m1 = R_PosInf, m2 = R_PosInf;
      for (int p = a0; p < a1; p++) {
        int s = nd->alive_idx[p];
        const double *rp = C_row + (size_t) s * K;
        double v = cost[s];
        for (int k = 0; k < K; k++) v += rp[k] * ctx->mu_best[k];
        if (v < m1) {
          m2 = m1;
          m1 = v;
        } else if (v < m2) {
          m2 = v;
        }
      }
      double cutoff = m1 + slack;
      for (int p = a0; p < a1; p++) {
        int s = nd->alive_idx[p];
        const double *rp = C_row + (size_t) s * K;
        double v = cost[s];
        for (int k = 0; k < K; k++) v += rp[k] * ctx->mu_best[k];
        if (v <= cutoff) {
          ctx->tmp_idx[new_total++] = s;
        }
      }
      if (new_total == ctx->tmp_off[h]) {
        free_node(nd);
        return NULL;
      }
      double g_h = m2 - m1;
      if (g_h <= slack && g_h < best_gap) {
        best_gap = g_h;
        branch_h = h;
      }
    }
    ctx->tmp_off[H] = new_total;

    /* Update node alive arrays in place (always shrinks: new_total <= total_alive) */
    memcpy(nd->alive_off, ctx->tmp_off, (size_t) (H + 1) * sizeof(int));
    if (new_total > 0) {
      memcpy(nd->alive_idx, ctx->tmp_idx, (size_t) new_total * sizeof(int));
    }
    nd->total_alive = new_total;
    nd->L = best_L;
    memcpy(nd->mu, ctx->mu_best, (size_t) K * sizeof(double));

    if (branch_h < 0) {
      /* All free groups collapsed to 1 state: loop to fix & improve */
      continue;
    }
    nd->branch = branch_h;
    return nd;
  }
}

SEXP C_menu_tree(SEXP cost_s, SEXP C_s, SEXP ends_s, SEXP lambda_s,
                 SEXP start_s, SEXP grain_s, SEXP max_nodes_s,
                 SEXP max_open_s, SEXP gap_s, SEXP exact_s) {
  if (!Rf_isReal(cost_s)) Rf_error("'cost' must be double");
  if (!Rf_isReal(C_s) || !Rf_isMatrix(C_s)) Rf_error("'C' must be a double matrix");
  if (!Rf_isInteger(ends_s)) Rf_error("'ends' must be integer");
  if (!Rf_isReal(lambda_s)) Rf_error("'lambda' must be double");

  const int N_orig = (int) XLENGTH(cost_s);
  const int K = Rf_ncols(C_s);
  const int H = LENGTH(ends_s);
  if (Rf_nrows(C_s) != N_orig || LENGTH(lambda_s) != K) {
    Rf_error("menu dimensions do not agree");
  }

  const double *cost_orig = REAL(cost_s);
  const double *C_orig = REAL(C_s);
  const int *ends = INTEGER(ends_s);
  for (int h = 0; h < H; h++) {
    int lo = (h == 0) ? 0 : ends[h - 1];
    if (ends[h] <= lo || ends[h] > N_orig) {
      Rf_error("'ends' must increase within 'cost'");
    }
  }
  const double *lam_init = REAL(lambda_s);
  const double grain = Rf_asReal(grain_s);
  const double max_nodes_d = Rf_asReal(max_nodes_s);
  const int max_nodes = (max_nodes_d > 2e9) ? 2000000000 : (int) max_nodes_d;
  const double max_open_d = Rf_asReal(max_open_s);
  const int max_open = (max_open_d > 2e9) ? 2000000000 : (int) max_open_d;
  const double gap = Rf_asReal(gap_s);
  const int exact = Rf_asLogical(exact_s) == TRUE;

  /* Ceiling cost and initial incumbent */
  double ceiling_cost = 1.0;
  for (int h = 0; h < H; h++) {
    int lo = (h == 0) ? 0 : ends[h - 1];
    int hi = ends[h];
    double mx = cost_orig[lo];
    for (int i = lo + 1; i < hi; i++) {
      if (cost_orig[i] > mx) mx = cost_orig[i];
    }
    ceiling_cost += mx;
  }

  int has_start = !Rf_isNull(start_s);
  double best_init = R_PosInf;
  int *start_orig = NULL;
  if (has_start) {
    if (!Rf_isInteger(start_s) || LENGTH(start_s) != H) {
      Rf_error("'start' must be an integer vector of length H");
    }
    start_orig = INTEGER(start_s);
    for (int h = 0; h < H; h++) {
      int lo = (h == 0) ? 0 : ends[h - 1];
      if (start_orig[h] == NA_INTEGER || start_orig[h] <= lo ||
          start_orig[h] > ends[h]) {
        Rf_error("'start' must pick a row of each group");
      }
    }
    best_init = 0.0;
    for (int h = 0; h < H; h++) best_init += cost_orig[start_orig[h] - 1];
  }
  double ref_scale = fabs(R_FINITE(best_init) ? best_init : ceiling_cost);
  if (ref_scale < 1e-300) ref_scale = 1e-300;
  double tol = 1e-9 * ref_scale;
  double root_lim = R_FINITE(best_init)
    ? (best_init - (grain - tol > tol ? grain - tol : tol))
    : ceiling_cost;

  /* Root Lagrangian + reduced-cost filtering + dominance */
  double *r_orig = (double *) malloc((size_t) (N_orig > 0 ? N_orig : 1) * sizeof(double));
  double *rmin_grp = (double *) malloc((size_t) H * sizeof(double));
  int *orig_to_new = (int *) malloc((size_t) (N_orig > 0 ? N_orig : 1) * sizeof(int));
  if (!r_orig || !rmin_grp || !orig_to_new) {
    free(r_orig); free(rmin_grp); free(orig_to_new);
    Rf_error("out of memory in C_menu_tree");
  }
  for (int i = 0; i < N_orig; i++) orig_to_new[i] = -1;

  double L_root = 0.0;
  for (int k = 0; k < K; k++) L_root -= lam_init[k];
  for (int h = 0; h < H; h++) {
    int lo = (h == 0) ? 0 : ends[h - 1];
    int hi = ends[h];
    double m1 = R_PosInf;
    for (int i = lo; i < hi; i++) {
      double v = cost_orig[i];
      for (int k = 0; k < K; k++) v += C_orig[i + (size_t) k * N_orig] * lam_init[k];
      r_orig[i] = v;
      if (v < m1) m1 = v;
    }
    rmin_grp[h] = m1;
    L_root += m1;
  }

  /* Filter & dominance per group, build compacted arrays */
  double root_slack = root_lim - L_root;
  int *kept_orig = (int *) malloc((size_t) (N_orig > 0 ? N_orig : 1) * sizeof(int));
  int *grp_off = (int *) malloc((size_t) (H + 1) * sizeof(int));
  SortItem *sbuf = (SortItem *) malloc((size_t) (N_orig > 0 ? N_orig : 1) * sizeof(SortItem));
  if (!kept_orig || !grp_off || !sbuf) {
    free(r_orig); free(rmin_grp); free(orig_to_new);
    free(kept_orig); free(grp_off); free(sbuf);
    Rf_error("out of memory in C_menu_tree");
  }

  int n_kept = 0;
  int any_empty = 0;
  for (int h = 0; h < H; h++) {
    grp_off[h] = n_kept;
    int lo = (h == 0) ? 0 : ends[h - 1];
    int hi = ends[h];
    int cnt = 0;
    double cutoff = rmin_grp[h] + root_slack;
    for (int i = lo; i < hi; i++) {
      if (r_orig[i] <= cutoff) {
        double cs = 0.0;
        for (int k = 0; k < K; k++) cs += C_orig[i + (size_t) k * N_orig];
        sbuf[cnt].orig = i;
        sbuf[cnt].cost = cost_orig[i];
        sbuf[cnt].csum = cs;
        cnt++;
      }
    }
    if (cnt > 1) {
      qsort(sbuf, (size_t) cnt, sizeof(SortItem), cmp_sort_item);
    }
    int g_start = n_kept;
    for (int a = 0; a < cnt; a++) {
      int i = sbuf[a].orig;
      int dom = 0;
      for (int b = g_start; b < n_kept && !dom; b++) {
        int j = kept_orig[b];
        int k = 0;
        while (k < K && C_orig[j + (size_t) k * N_orig] <=
                        C_orig[i + (size_t) k * N_orig] + 1e-15) {
          k++;
        }
        if (k == K) dom = 1;
      }
      if (!dom) {
        orig_to_new[i] = n_kept;
        kept_orig[n_kept++] = i;
      }
    }
    if (n_kept == g_start) any_empty = 1;
  }
  grp_off[H] = n_kept;
  free(r_orig);
  free(rmin_grp);
  free(sbuf);

  /* Initialize MenuCtx */
  MenuCtx ctx;
  memset(&ctx, 0, sizeof(ctx));
  ctx.H = H;
  ctx.K = K;
  ctx.N = n_kept;
  ctx.kept_total = n_kept;
  ctx.grain = grain;
  ctx.tol = tol;
  ctx.ceiling_cost = ceiling_cost;
  ctx.best = best_init;
  ctx.has_plan = 0;
  ctx.exact = exact;
  ctx.grp_off = grp_off;

  ctx.cost = (double *) malloc((size_t) (n_kept > 0 ? n_kept : 1) * sizeof(double));
  ctx.C_row = (double *) malloc((size_t) (n_kept > 0 ? n_kept : 1) * (size_t) (K > 0 ? K : 1) * sizeof(double));
  ctx.orig_idx = (int *) malloc((size_t) (n_kept > 0 ? n_kept : 1) * sizeof(int));
  ctx.best_plan = (int *) malloc((size_t) H * sizeof(int));
  ctx.free_list = (int *) malloc((size_t) H * sizeof(int));
  ctx.cmin = (double *) malloc((size_t) (K > 0 ? K : 1) * sizeof(double));
  ctx.mu_cur = (double *) malloc((size_t) (K > 0 ? K : 1) * sizeof(double));
  ctx.mu_best = (double *) malloc((size_t) (K > 0 ? K : 1) * sizeof(double));
  ctx.grad = (double *) malloc((size_t) (K > 0 ? K : 1) * sizeof(double));
  ctx.grp_arg = (int *) malloc((size_t) H * sizeof(int));
  ctx.tmp_off = (int *) malloc((size_t) (H + 1) * sizeof(int));
  ctx.tmp_idx = (int *) malloc((size_t) (n_kept > 0 ? n_kept : 1) * sizeof(int));
  ctx.tmp_plan = (int *) malloc((size_t) H * sizeof(int));

  if (!ctx.cost || !ctx.C_row || !ctx.orig_idx || !ctx.best_plan ||
      !ctx.free_list || !ctx.cmin || !ctx.mu_cur || !ctx.mu_best ||
      !ctx.grad || !ctx.grp_arg ||
      !ctx.tmp_off || !ctx.tmp_idx || !ctx.tmp_plan) {
    free(kept_orig);
    free(orig_to_new);
    free_ctx(&ctx);
    Rf_error("out of memory in C_menu_tree");
  }

  for (int s = 0; s < n_kept; s++) {
    int orig = kept_orig[s];
    ctx.cost[s] = cost_orig[orig];
    ctx.orig_idx[s] = orig + 1; /* 1-based */
    for (int k = 0; k < K; k++) {
      ctx.C_row[(size_t) s * K + k] = C_orig[orig + (size_t) k * N_orig];
    }
  }
  free(kept_orig);
  free(orig_to_new);

  /* Build and evaluate root */
  int nodes = 1;
  int interrupted = 0;
  if (!any_empty) {
    MenuNode *root = alloc_node(H, K, n_kept);
    if (!root) {
      free_ctx(&ctx);
      Rf_error("out of memory in C_menu_tree");
    }
    for (int h = 0; h < H; h++) root->plan[h] = -1;
    for (int k = 0; k < K; k++) {
      root->cap[k] = 1.0;
      root->mu[k] = lam_init[k];
    }
    memcpy(root->alive_off, grp_off, (size_t) (H + 1) * sizeof(int));
    for (int s = 0; s < n_kept; s++) root->alive_idx[s] = s;
    root->L = L_root;
    root->c_fix = 0.0;

    root = evaluate_node(&ctx, root);
    if (root) {
      if (!heap_push(&ctx, root)) {
        free_node(root);
        free_ctx(&ctx);
        Rf_error("out of memory in C_menu_tree");
      }
    }
  }

  /* Branch-and-bound main loop */
  int iter = 0;
  while (ctx.open_size > 0 || ctx.stack_size > 0) {
    if (nodes >= max_nodes) break;
    if ((++iter & 255) == 0 && check_interrupted()) {
      interrupted = 1;
      break;
    }
    if (gap > 0.0 && R_FINITE(ctx.best)) {
      double ob = ctx.best;
      if (ctx.open_size > 0 && ctx.open[0]->L < ob) ob = ctx.open[0]->L;
      for (int i = 0; i < ctx.stack_size; i++) {
        if (ctx.stack[i]->L < ob) ob = ctx.stack[i]->L;
      }
      if (ctx.best - ob <= gap * fabs(ctx.best)) break;
    }

    MenuNode *nd = NULL;
    int deep = 0;
    int from_stack = 0;
    if (ctx.stack_size > 0) {
      nd = ctx.stack[--ctx.stack_size];
      if (nd->L > ctx_limit(&ctx)) {
        free_node(nd);
        continue;
      }
      deep = 1;
      from_stack = 1;
    } else {
      if (ctx.open[0]->L > ctx_limit(&ctx)) {
        for (int i = 0; i < ctx.open_size; i++) free_node(ctx.open[i]);
        ctx.open_size = 0;
        break;
      }
      nd = heap_pop(&ctx);
      deep = (ctx.open_size >= max_open);
    }

    int bh = nd->branch;
    int b0 = nd->alive_off[bh];
    int b1 = nd->alive_off[bh + 1];
    int n_branch = b1 - b0;
    if ((double) nodes + n_branch > max_nodes) {
      int ok = from_stack ? stack_push(&ctx, nd) : heap_push(&ctx, nd);
      if (!ok) {
        free_node(nd);
        free_ctx(&ctx);
        Rf_error("out of memory in C_menu_tree");
      }
      break;
    }
    if (n_branch > ctx.kids_cap) {
      MenuNode **tmp = (MenuNode **) realloc(ctx.kids, (size_t) n_branch * sizeof(MenuNode *));
      if (!tmp) {
        free_node(nd);
        free_ctx(&ctx);
        Rf_error("out of memory in C_menu_tree");
      }
      ctx.kids = tmp;
      ctx.kids_cap = n_branch;
    }

    int n_kids = 0;
    int child_alive_len = nd->total_alive - n_branch + 1;
    for (int bi = b0; bi < b1; bi++) {
      int s_pick = nd->alive_idx[bi];
      MenuNode *ch = alloc_node(H, K, child_alive_len);
      if (!ch) {
        for (int j = 0; j < n_kids; j++) free_node(ctx.kids[j]);
        free_node(nd);
        free_ctx(&ctx);
        Rf_error("out of memory in C_menu_tree");
      }
      memcpy(ch->plan, nd->plan, (size_t) H * sizeof(int));
      memcpy(ch->cap, nd->cap, (size_t) K * sizeof(double));
      memcpy(ch->mu, nd->mu, (size_t) K * sizeof(double));
      ch->c_fix = nd->c_fix;
      ch->L = nd->L;

      int pos = 0;
      for (int h = 0; h < H; h++) {
        ch->alive_off[h] = pos;
        if (h == bh) {
          ch->alive_idx[pos++] = s_pick;
        } else {
          int len = nd->alive_off[h + 1] - nd->alive_off[h];
          if (len > 0) {
            memcpy(ch->alive_idx + pos, nd->alive_idx + nd->alive_off[h],
                   (size_t) len * sizeof(int));
            pos += len;
          }
        }
      }
      ch->alive_off[H] = pos;
      ch->total_alive = pos;

      nodes++;
      ch = evaluate_node(&ctx, ch);
      if (ch) {
        if (ch->L <= ctx_limit(&ctx)) {
          ctx.kids[n_kids++] = ch;
        } else {
          free_node(ch);
        }
      }
    }
    free_node(nd);

    if (n_kids == 0) continue;
    if (deep) {
      if (n_kids > 1) {
        sort_kids_desc(ctx.kids, n_kids);
      }
      for (int j = 0; j < n_kids; j++) {
        if (!stack_push(&ctx, ctx.kids[j])) {
          for (int r = j; r < n_kids; r++) free_node(ctx.kids[r]);
          free_ctx(&ctx);
          Rf_error("out of memory in C_menu_tree");
        }
      }
    } else {
      for (int j = 0; j < n_kids; j++) {
        if (!heap_push(&ctx, ctx.kids[j])) {
          for (int r = j; r < n_kids; r++) free_node(ctx.kids[r]);
          free_ctx(&ctx);
          Rf_error("out of memory in C_menu_tree");
        }
      }
    }
  }

  if (interrupted) {
    free_ctx(&ctx);
    Rf_error("the whole-unit search was interrupted");
  }

  /* Summarize remaining open/stack bounds */
  double lim_final = ctx_limit(&ctx);
  int n_left = 0;
  double min_left = R_PosInf;
  for (int i = 0; i < ctx.open_size; i++) {
    double lv = ctx.open[i]->L;
    if (lv <= lim_final) {
      n_left++;
      if (lv < min_left) min_left = lv;
    }
  }
  for (int i = 0; i < ctx.stack_size; i++) {
    double lv = ctx.stack[i]->L;
    if (lv <= lim_final) {
      n_left++;
      if (lv < min_left) min_left = lv;
    }
  }

  SEXP pick_sexp = R_NilValue;
  double out_val = R_PosInf;
  double out_bound = R_PosInf;
  const char *status_str = "infeasible";

  if (!R_FINITE(ctx.best)) {
    if (n_left == 0) {
      status_str = "infeasible";
      out_bound = R_PosInf;
    } else {
      status_str = "limit";
      out_bound = min_left;
    }
  } else {
    pick_sexp = PROTECT(Rf_allocVector(INTSXP, H));
    int *pp = INTEGER(pick_sexp);
    if (ctx.has_plan) {
      for (int h = 0; h < H; h++) pp[h] = ctx.orig_idx[ctx.best_plan[h]];
    } else {
      for (int h = 0; h < H; h++) pp[h] = start_orig[h];
    }
    out_val = ctx.best;
    out_bound = (n_left > 0) ? (min_left < ctx.best ? min_left : ctx.best) : ctx.best;
    double g_thr = gap > 1e-6 ? gap : 1e-6;
    int near = (ctx.best - out_bound <= g_thr * fabs(ctx.best));
    status_str = (n_left == 0) ? "optimal" : (near ? "gap" : "limit");
  }

  int n_prot = R_FINITE(ctx.best) ? 1 : 0;
  int kept_final = ctx.kept_total;
  free_ctx(&ctx);

  SEXP out = PROTECT(Rf_allocVector(VECSXP, 6));
  SET_VECTOR_ELT(out, 0, pick_sexp);
  SET_VECTOR_ELT(out, 1, Rf_ScalarReal(out_val));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(out_bound));
  SET_VECTOR_ELT(out, 3, Rf_mkString(status_str));
  SET_VECTOR_ELT(out, 4, Rf_ScalarReal((double) nodes));
  SET_VECTOR_ELT(out, 5, Rf_ScalarInteger(kept_final));

  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 6));
  SET_STRING_ELT(nms, 0, Rf_mkChar("pick"));
  SET_STRING_ELT(nms, 1, Rf_mkChar("value"));
  SET_STRING_ELT(nms, 2, Rf_mkChar("bound"));
  SET_STRING_ELT(nms, 3, Rf_mkChar("status"));
  SET_STRING_ELT(nms, 4, Rf_mkChar("nodes"));
  SET_STRING_ELT(nms, 5, Rf_mkChar("kept"));
  Rf_setAttrib(out, R_NamesSymbol, nms);

  UNPROTECT(n_prot + 2);
  return out;
}
