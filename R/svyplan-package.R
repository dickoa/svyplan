#' @keywords internal
#'
#' @section Where svyplan stops:
#' svyplan stops at the plan. It decides how many units, allocated where, to
#' what precision, and it never draws a sample: no selection probabilities,
#' no weights, no drawn units come out of any function here. What you get is
#' the design a sampler is then asked to realize.
#'
#' When the plan is settled, \pkg{sondage} draws it. Analysis of the realized
#' sample belongs to \pkg{survey} or \pkg{srvyr}.
#'
#' @section Nonresponse adjustment:
#' Dividing by `resp_rate` is an expected-information calculation. It sets the
#' expected respondent count and evaluates variance at that net size. It
#' assumes response is ignorable under the adjustment planned for analysis.
#' No sample-size inflation removes nonresponse bias, and this calculation
#' does not include variance from response weights. Informative nonresponse
#' calls for modeling, weighting, follow-up design, or sensitivity analysis
#' over `resp_rate`. [n_twophase()] provides the package's explicit
#' nonresponse follow-up design, and [predict()] supports sensitivity grids.
#'
#' @section Which function:
#' Eleven problems, each with the function or pair that answers it. Every
#' `n_` has a `prec_` reading the same design back the other way, so a row
#' names one thing to learn rather than two.
#'
#' \tabular{ll}{
#'   **I want to size, or evaluate** \tab **Functions** \cr
#'   a proportion, a mean or a ratio of two totals \tab `n_prop()`/`prec_prop()`, `n_mean()`/`prec_mean()`, `n_ratio()`/`prec_ratio()` \cr
#'   several indicators at once, taking the most demanding \tab `n_multi()`/`prec_multi()` \cr
#'   a two- or three-stage cluster design, for one indicator or a table of them \tab `n_cluster()`/`prec_cluster()` \cr
#'   an allocation across strata or domains, under a budget or a CV target \tab `n_alloc()`/`prec_alloc()`, `strata_bound()` \cr
#'   a two-phase design that screens or follows up \tab `n_twophase()`/`prec_twophase()` \cr
#'   a change between two occasions of a repeated survey \tab `n_change()`/`prec_change()` \cr
#'   the average of several occasions of a repeated survey \tab `n_pooled()`/`prec_pooled()` \cr
#'   a panel that must still deliver a sample after attrition \tab `n_panel()`/`prec_panel()` \cr
#'   the overlap and field schedule a rotation pattern produces \tab `design_rotation()`, `design_overlap()`, `design_schedule()` \cr
#'   the power of a two-group comparison or a difference-in-differences \tab `power_prop()`, `power_mean()`, `power_did()` \cr
#'   the design effect, effective size or degrees of freedom of a plan \tab `design_effect()`, `effective_n()`, `design_df()`, `varcomp()`
#' }
#'
#' @section Notation:
#' Argument names spell out what a quantity is rather than reproducing the
#' symbol used in any one textbook. Readers coming from the standard
#' references can map them as follows.
#'
#' \tabular{lll}{
#'   **Argument** \tab **Symbol** \tab **Meaning** \cr
#'   `icc` \tab \eqn{\delta} \tab Design-based measure of homogeneity within
#'     clusters, \eqn{V_b/(V_b+V_w)}. Constrained to \eqn{[0, 1]}, so it is
#'     not interchangeable with a mixed-model ICC, which can be negative.
#'     Written \eqn{\delta} by Valliant, Dever and Kreuter (2018) and
#'     related to Kish's \emph{roh}. \cr
#'   `var_ratio` \tab \eqn{k} \tab Ratio of the stage components' unit
#'     variance to the analysis variable's. Defaults to 1. \cr
#'   `unit_relvar` \tab \eqn{V} \tab Unit relvariance, \eqn{S^2/\bar{y}^2}{S^2/ybar^2},
#'     that is the squared population coefficient of variation. \cr
#'   `deff` \tab \eqn{DEFF} \tab Design effect. \cr
#'   `n_psu` \tab \eqn{n_1} \tab Number of PSUs selected. \cr
#'   `n_per_psu` \tab \eqn{n_2} \tab Units selected \emph{per} PSU, which is
#'     a sample take, not the PSU's population size. \cr
#'   `n_per_ssu` \tab \eqn{n_3} \tab Units selected per SSU. \cr
#'   `moe` \tab \eqn{e} \tab Margin of error, the half-width of the
#'     confidence interval. \cr
#'   `cv` \tab \eqn{CV} \tab Coefficient of variation, the relative standard
#'     error. \cr
#'   `overlap`, `overlap_cor` \tab \eqn{\gamma}, \eqn{\rho} \tab Panel
#'     overlap fraction and the correlation between occasions. \cr
#'   `alloc_q` \tab \eqn{q} \tab Bankier power-allocation exponent, used only
#'     when `alloc = "power"`. It is unrelated to statistical power. \cr
#' }
#'
#' Dispersion may be given as either `var` or `sd` wherever a mean is being
#' planned, and exactly one is required. The exception is [n_twophase()],
#' whose frame
#' takes `sd` only.
#'
#' @references
#' Cochran, W. G. (1977). \emph{Sampling Techniques}, 3rd edition. Wiley.
#'
#' Kish, L. (1965). \emph{Survey Sampling}. Wiley.
#'
#' Valliant, R., Dever, J. A., and Kreuter, F. (2018). \emph{Practical Tools
#' for Designing and Weighting Survey Samples}, 2nd edition. Springer.
"_PACKAGE"

#' @importFrom graphics abline axis barplot hist mtext par plot.new
#'   plot.window points rect strwidth text title
#' @importFrom grDevices adjustcolor nclass.FD
#' @importFrom stats lm.fit optim optimize plogis pnorm predict qlogis quantile qnorm sd terms uniroot var weighted.mean setNames
#' @importFrom utils modifyList
NULL
