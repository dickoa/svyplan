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
