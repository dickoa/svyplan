/* One part of a stratum's menu of field designs, listed by an event sweep.
 *
 * A part fixes the PSUs forced into certainty. Over a range of sampling
 * fractions f it lists every distinct state the closure gives: the
 * threshold test of .psu_classify(), the whole-PSU draw of
 * .psu_stratum_counts() and the crossing test of .psu_crossing(), zone by
 * zone when zones are drawn. Each is evaluated with the float expressions
 * of the R code, in its order, and with its sums accumulated as R's sum()
 * and cumsum() accumulate them, so a listed state is the state the R code
 * computes at that f.
 *
 * The state is constant between consecutive events. The sweep lists the
 * state at f_lo, then for each event f the open interval to the next
 * event, at its midpoint, and the next event itself. The next event is
 * read from the closure just right of f. */
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

/* Rounds a product before it is added to, so no compiler fuses the two
 * into one multiply-add. Fusing changes the last bit on some platforms. */
static double rounded(double x) {
  volatile double v = x;
  return v;
}

typedef struct part part_t;

struct part {
  int n;               /* PSUs */
  const double *M;     /* sizes */
  double N;            /* stratum size, from the frame */
  double b;            /* take per remainder PSU */
  double cutoff;
  const int *forced;   /* logical */
  int m;               /* PSUs per zone, 0 without zones */
  const int *ord;      /* zone order of every PSU, 0-based */
  int *enter;          /* round each PSU entered certainty, INT_MAX if not */
  int K;               /* rounds of the last closure */
  double last_Q;       /* certainty take and remainder draw of that closure */
  double last_a;
  double *round_Q;     /* certainty takes of every closure round */
  int *round_avail;    /* remainder counts of every closure round */
  int *hit;
  int *zone;
  void *zacc;          /* zone totals in the accumulator's type */
  double *zsum;
  int *zcount;
  int *cur;            /* scratch for the final set */
  /* The copies below for the accumulator R's sum() and cumsum() would
   * give the same result with: double when every size is whole and their
   * total below 2^53, where both sum exactly, long double otherwise where
   * R is built with it. */
  double (*draw_of)(const part_t *, int, double, double *, int *);
  int (*crossing)(part_t *, int, double);
};

static int in_set(const part_t *p, int i, int r) {
  return p->enter[i] <= r;
}

#define ACC_T double
#define SUFFIX(name) name##_double
#include "psu_menu_sums.h"
#undef ACC_T
#undef SUFFIX

#define ACC_T long double
#define SUFFIX(name) name##_long
#include "psu_menu_sums.h"
#undef ACC_T
#undef SUFFIX

/* The closure at f, recorded as the round each PSU entered. */
static void closure(part_t *p, double f) {
  double n_h = f * p->N;
  double thr = p->cutoff * p->b / (n_h / p->N);
  for (int i = 0; i < p->n; i++) {
    p->enter[i] = (p->M[i] >= thr || p->forced[i]) ? 0 : INT_MAX;
  }
  p->K = 0;
  for (;;) {
    double k = p->draw_of(p, p->K, n_h, &p->last_Q,
                          &p->round_avail[p->K]);
    p->last_a = k;
    p->round_Q[p->K] = p->last_Q;
    if (!p->crossing(p, p->K, k)) break;
    p->K++;
    for (int i = 0; i < p->n; i++) {
      if (p->hit[i]) p->enter[i] = p->K;
    }
  }
}

/* The next event above f of the closure read at fr, just right of f: a take
 * j = ceil(fr M_i) changes just after j / M_i, a draw k = ceil((fr N - Q)
 * / b) just after (k b + Q) / N, and C_0 at a threshold cutoff b / M_i.
 * No interval longer than fr - f is skipped. */
static double next_event(part_t *p, double f, double fr) {
  double e = R_PosInf;
  double n_h = fr * p->N;
  for (int i = 0; i < p->n; i++) {
    if (in_set(p, i, p->K)) {
      double j = ceil(n_h / p->N * p->M[i]);
      if (j < p->M[i]) {
        double t = j / p->M[i];
        if (t > f && t < e) e = t;
      }
    }
  }
  for (int r = 0; r <= p->K; r++) {
    if (p->round_avail[r] == 0) continue;
    double nc = p->round_Q[r];
    double rest = n_h - nc;
    if (rest < 0) rest = 0;
    double k = ceil(rest / p->b);
    double t = (rounded(k * p->b) + nc) / p->N;
    if (t > f && t < e) e = t;
  }
  for (int i = 0; i < p->n; i++) {
    if (p->enter[i] != 0 && !p->forced[i]) {
      double t = p->cutoff * p->b / p->M[i];
      if (t > f && t < e) e = t;
    }
  }
  return e;
}

/* Output columns, grown in R_alloc memory so an error or an interrupt
 * leaks nothing. */
typedef struct {
  R_xlen_t len, cap;
  double *f, *Q, *A, *held, *n_c, *a;
  int *set;
  int n_sets, cap_sets;
  int *bits;           /* n x cap_sets */
} out_t;

/* R_alloc() memory is aligned for double only, and long double needs 16
 * bytes on x86_64, so its arrays are over-allocated and aligned here. */
static void *aligned_alloc_r(size_t count, size_t size, size_t align) {
  char *raw = R_alloc(count * size + align, 1);
  uintptr_t at = (uintptr_t) raw;
  at = (at + align - 1) & ~(uintptr_t) (align - 1);
  return (void *) at;
}

static void *grow(void *old, size_t used, size_t cap, size_t size) {
  void *v = R_alloc(cap, (int) size);
  if (used) memcpy(v, old, used * size);
  return v;
}

static void out_reserve(out_t *o) {
  if (o->len < o->cap) return;
  R_xlen_t cap = o->cap ? 2 * o->cap : 1024;
  o->f = grow(o->f, o->len, cap, sizeof(double));
  o->Q = grow(o->Q, o->len, cap, sizeof(double));
  o->A = grow(o->A, o->len, cap, sizeof(double));
  o->held = grow(o->held, o->len, cap, sizeof(double));
  o->n_c = grow(o->n_c, o->len, cap, sizeof(double));
  o->a = grow(o->a, o->len, cap, sizeof(double));
  o->set = grow(o->set, o->len, cap, sizeof(int));
  o->cap = cap;
}

/* Index of the final set of the last closure, added when new. */
static int set_id(const part_t *p, out_t *o) {
  int *cur = p->cur;
  for (int i = 0; i < p->n; i++) cur[i] = in_set(p, i, p->K);
  for (int k = o->n_sets - 1; k >= 0; k--) {
    if (memcmp(o->bits + (size_t) k * p->n, cur, sizeof(int) * p->n) == 0) return k;
  }
  if (o->n_sets == o->cap_sets) {
    int cap = o->cap_sets ? 2 * o->cap_sets : 16;
    o->bits = grow(o->bits, (size_t) o->n_sets * p->n, (size_t) cap * p->n, sizeof(int));
    o->cap_sets = cap;
  }
  memcpy(o->bits + (size_t) o->n_sets * p->n, cur, sizeof(int) * p->n);
  return o->n_sets++;
}

/* The state at f: the sum of its takes, as the draw sums them, sum(M^2 / q)
 * and the total size over its certainty PSUs, and its draw. */
static void push(part_t *p, double f, out_t *o) {
  closure(p, f);
  double n_h = f * p->N;
  double Q = p->last_Q;
  double a = p->last_a;
  int set = set_id(p, o) + 1;
  if (o->len > 0) {
    R_xlen_t last = o->len - 1;
    if (o->set[last] == set && o->Q[last] == Q && o->a[last] == a) return;
  }
  double A = 0, held = 0;
  int n_c = 0;
  for (int i = 0; i < p->n; i++) {
    if (in_set(p, i, p->K)) {
      double q = ceil(n_h / p->N * p->M[i]);
      if (q > p->M[i]) q = p->M[i];
      A += p->M[i] * p->M[i] / q;
      held += p->M[i];
      n_c++;
    }
  }
  out_reserve(o);
  R_xlen_t k = o->len++;
  o->f[k] = f;
  o->Q[k] = Q;
  o->A[k] = A;
  o->held[k] = held;
  o->n_c[k] = n_c;
  o->a[k] = a;
  o->set[k] = set;
}

static void check_real(SEXP x, R_xlen_t len, const char *what) {
  if (TYPEOF(x) != REALSXP || (len >= 0 && XLENGTH(x) != len)) {
    Rf_error("internal error: '%s' must be a double vector of length %ld",
             what, (long) len);
  }
}

SEXP C_psu_menu_part(SEXP sizes, SEXP N, SEXP take, SEXP f_lo, SEXP f_hi,
                     SEXP cutoff, SEXP forced, SEXP m, SEXP ord,
                     SEXP f_stop, SEXP ld) {
  R_xlen_t n = XLENGTH(sizes);
  if (n < 1 || n > INT_MAX / 4) Rf_error("internal error: a part needs 1 to %d PSUs", INT_MAX / 4);
  check_real(sizes, n, "sizes");
  check_real(N, 1, "N");
  check_real(take, 1, "take");
  check_real(f_lo, 1, "f_lo");
  check_real(f_hi, 1, "f_hi");
  check_real(cutoff, 1, "cutoff");
  check_real(f_stop, 1, "f_stop");
  if (TYPEOF(forced) != LGLSXP || XLENGTH(forced) != n) {
    Rf_error("internal error: 'forced' must be a logical vector, one per PSU");
  }
  if (TYPEOF(ord) != INTSXP || XLENGTH(ord) != n) {
    Rf_error("internal error: 'ord' must be an integer vector, one per PSU");
  }
  if (TYPEOF(m) != INTSXP || XLENGTH(m) != 1 || TYPEOF(ld) != LGLSXP || XLENGTH(ld) != 1) {
    Rf_error("internal error: 'm' must be one integer and 'ld' one logical");
  }

  part_t p;
  p.n = (int) n;
  p.M = REAL(sizes);
  p.N = REAL(N)[0];
  p.b = REAL(take)[0];
  p.cutoff = REAL(cutoff)[0];
  p.forced = LOGICAL(forced);
  p.m = INTEGER(m)[0];
  p.ord = INTEGER(ord);
  int whole = 1;
  double total = 0;
  for (int i = 0; i < p.n; i++) {
    if (!R_FINITE(p.M[i]) || p.M[i] <= 0 || p.forced[i] == NA_LOGICAL ||
        p.ord[i] < 0 || p.ord[i] >= p.n) {
      Rf_error("internal error: invalid PSU %d", i + 1);
    }
    if (p.M[i] != floor(p.M[i])) whole = 0;
    total += p.M[i];
  }
  if (whole && total < 9007199254740992.0) {
    p.draw_of = draw_of_double;
    p.crossing = crossing_double;
  } else if (LOGICAL(ld)[0] == TRUE) {
    p.draw_of = draw_of_long;
    p.crossing = crossing_long;
  } else {
    p.draw_of = draw_of_double;
    p.crossing = crossing_double;
  }
  double lo = REAL(f_lo)[0], hi = REAL(f_hi)[0], stop_at = REAL(f_stop)[0];
  if (!(R_FINITE(p.N) && p.N > 0 && R_FINITE(p.b) && p.b > 0 &&
        p.cutoff > 0 && p.cutoff <= 1 && p.m >= 0 && R_FINITE(lo) &&
        lo > 0 && R_FINITE(hi) && hi >= lo && hi <= 1 && !ISNAN(stop_at))) {
    Rf_error("internal error: invalid part settings");
  }
  p.enter = (int *) R_alloc(p.n, sizeof(int));
  p.hit = (int *) R_alloc(p.n, sizeof(int));
  p.zone = (int *) R_alloc(p.n, sizeof(int));
  p.zacc = aligned_alloc_r((size_t) p.n, sizeof(long double), _Alignof(long double));
  p.zsum = (double *) R_alloc(p.n, sizeof(double));
  p.cur = (int *) R_alloc(p.n, sizeof(int));
  p.round_Q = (double *) R_alloc((size_t) (p.n + 1), sizeof(double));
  p.round_avail = (int *) R_alloc((size_t) (p.n + 1), sizeof(int));
  p.zcount = (int *) R_alloc(p.n, sizeof(int));

  out_t o;
  memset(&o, 0, sizeof(o));
  double f = lo;
  int stopped = 0;
  unsigned long steps = 0;
  push(&p, f, &o);
  while (f < hi) {
    if (f >= stop_at) {
      stopped = 1;
      break;
    }
    if ((++steps & 4095UL) == 0) R_CheckUserInterrupt();
    /* Just right of f, past the few ulps the closure's float expressions
       err by. An event closer to f than that is merged with it. */
    double fr = f * (1 + 64 * DBL_EPSILON);
    closure(&p, fr);
    double e = next_event(&p, f, fr);
    if (!(e > f)) Rf_error("internal error: the menu sweep did not advance at f = %g", f);
    double mid = R_FINITE(e) ? 0.5 * (f + e) : fr;
    /* The state is constant up to e, so a point inside the part stands for
       it. Past hi the fraction would misreport the state's threshold. */
    if (mid > hi) mid = 0.5 * (f + hi);
    if (mid > f && (!R_FINITE(e) || mid < e)) push(&p, mid, &o);
    if (!R_FINITE(e) || e > hi) break;
    push(&p, e, &o);
    f = e;
  }

  const char *names[] = {"f", "set", "Q", "A", "held", "n_c", "a", ""};
  SEXP states = PROTECT(Rf_mkNamed(VECSXP, names));
  double *cols[] = {o.f, NULL, o.Q, o.A, o.held, o.n_c, o.a};
  for (int j = 0; j < 7; j++) {
    SEXP v;
    if (j == 1) {
      v = PROTECT(Rf_allocVector(INTSXP, o.len));
      if (o.len) memcpy(INTEGER(v), o.set, sizeof(int) * o.len);
    } else {
      v = PROTECT(Rf_allocVector(REALSXP, o.len));
      if (o.len) memcpy(REAL(v), cols[j], sizeof(double) * o.len);
    }
    SET_VECTOR_ELT(states, j, v);
    UNPROTECT(1);
  }
  SEXP sets = PROTECT(Rf_allocMatrix(LGLSXP, p.n, o.n_sets));
  int *sl = LOGICAL(sets);
  for (size_t k = 0; k < (size_t) o.n_sets * p.n; k++) sl[k] = o.bits[k];
  const char *out_names[] = {"states", "sets", "f_end", "stopped", ""};
  SEXP out = PROTECT(Rf_mkNamed(VECSXP, out_names));
  SET_VECTOR_ELT(out, 0, states);
  SET_VECTOR_ELT(out, 1, sets);
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(f));
  SET_VECTOR_ELT(out, 3, Rf_ScalarLogical(stopped));
  UNPROTECT(3);
  return out;
}
