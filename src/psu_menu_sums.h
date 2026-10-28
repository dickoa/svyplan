/* The functions of psu_menu.c that accumulate sums, written once for each
 * accumulator. psu_menu.c includes this file twice, with ACC_T double and
 * long double and SUFFIX naming each copy. */

/* The draw of round r's set at n_h units, .psu_stratum_counts(). */
static double SUFFIX(draw_of)(const part_t *p, int r, double n_h,
                              double *n_c_int, int *avail_out) {
  ACC_T nc = 0;
  int avail = 0;
  for (int i = 0; i < p->n; i++) {
    if (in_set(p, i, r)) {
      double q = ceil(n_h / p->N * p->M[i]);
      nc += q < p->M[i] ? q : p->M[i];
    } else {
      avail++;
    }
  }
  double nci = (double) nc;
  if (n_c_int) *n_c_int = nci;
  if (avail_out) *avail_out = avail;
  double rest = n_h - nci;
  if (rest < 0) rest = 0;
  double k = ceil(rest / p->b);
  double floor_k = avail > 0 ? 1 : 0;
  if (k < floor_k) k = floor_k;
  if (p->m > 1 && k > 0) k = p->m * ceil(k / p->m);
  if (k > avail) k = avail;
  return k;
}

/* PSUs of round r's remainder the draw k takes whole, .psu_crossing():
 * the stratum-wide test, then, when it takes none and zones are drawn,
 * the zones cut on cumulative size in the zone order, .psu_zones(). */
static int SUFFIX(crossing)(part_t *p, int r, double k) {
  int any = 0;
  ACC_T rs = 0;
  int n_rest = 0;
  for (int i = 0; i < p->n; i++) {
    p->hit[i] = 0;
    if (!in_set(p, i, r)) {
      rs += p->M[i];
      n_rest++;
    }
  }
  if (k <= 0 || n_rest == 0) return 0;
  double rest_size = (double) rs;
  double level = p->cutoff - 100 * DBL_EPSILON;
  for (int i = 0; i < p->n; i++) {
    if (!in_set(p, i, r) && k * p->M[i] / rest_size >= level) {
      p->hit[i] = 1;
      any = 1;
    }
  }
  if (p->m < 1 || any) return any;
  double nz = floor(k / p->m);
  if (nz < 1) return 0;
  int n_zone = (int) nz;
  /* The last cumulative size is the total the width divides. */
  ACC_T cum = 0;
  for (int j = 0; j < p->n; j++) {
    int i = p->ord[j];
    if (!in_set(p, i, r)) cum += p->M[i];
  }
  double width = (double) cum / n_zone;
  cum = 0;
  for (int j = 0; j < p->n; j++) {
    int i = p->ord[j];
    if (in_set(p, i, r)) continue;
    cum += p->M[i];
    double z = ceil((double) cum / width - 1e-9);
    if (z < 1) z = 1;
    if (z > n_zone) z = n_zone;
    p->zone[i] = (int) z - 1;
  }
  /* Zone totals in index order, as sum(sizes[which(zone == z)]), in one
   * pass: each zone's terms still come in increasing index. */
  ACC_T *zacc = (ACC_T *) p->zacc;
  for (int z = 0; z < n_zone; z++) {
    zacc[z] = 0;
    p->zcount[z] = 0;
  }
  for (int i = 0; i < p->n; i++) {
    if (in_set(p, i, r)) continue;
    zacc[p->zone[i]] += p->M[i];
    p->zcount[p->zone[i]]++;
  }
  for (int z = 0; z < n_zone; z++) p->zsum[z] = (double) zacc[z];
  for (int i = 0; i < p->n; i++) {
    if (in_set(p, i, r)) continue;
    int z = p->zone[i];
    if (p->zcount[z] <= p->m || p->m * p->M[i] / p->zsum[z] >= level) {
      p->hit[i] = 1;
      any = 1;
    }
  }
  return any;
}
