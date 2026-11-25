# svyplan 0.26.2

Initial release of the package.

- Size or assess proportions, means, and ratios using precision targets,
  finite populations, design effects, eligibility, and response assumptions.
  Size several indicators or domains together and check respondent floors.
- Plan stratified and multistage cluster designs. Allocate a fixed total,
  seek a precision target, or work within a budget. Cluster designs support
  stage-specific costs, response rates, and finite population corrections.
- Use `n_alloc_joint()` and `prec_alloc_joint()` for joint indicator and
  domain requirements, budget objectives that minimize a weighted sum of
  relative variances or the largest of them, and PSU registers with certainty
  units and optional zones. Read the recommended whole design and search
  status alongside the continuous calculation.
- Plan two-phase sampling and nonresponse follow-up, compare repeated
  estimates, and size panel recruitment. Describe rotations, overlap, and
  fieldwork schedules.
- Calculate power, construct strata, estimate variance components, and
  obtain design effects, effective sample sizes, and degrees of freedom.
- Reuse assumptions in `svyplan()` profiles, assess fits with paired
  precision functions, and compare scenarios using `predict()`.
