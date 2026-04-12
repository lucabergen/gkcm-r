
#' GKCM conditional independence test
#'
#' Tests whether variable `x` is conditionally independent of variable `y` given
#' conditioning variables `Z` using the Generalized Kernel Covariance Measure (GKCM).
#' Supports mixed continuous and categorical data and returns a p-value together with
#' a standardized dependence measure.
#'
#' @details
#' `x` and `y` are embedded using a Gaussian kernel for numeric variables (bandwidth
#' \eqn{\sigma = 1}) and a Dirac kernel for categorical variables.
#'
#' Set `Z = NULL` to perform a marginal (unconditional) independence test. If a
#' 0-column `data.frame` or matrix is supplied for `Z`, it is interpreted as `Z = NULL`
#' and a warning is issued.
#'
#' For causal discovery workflows (e.g., PC/FCI), use [gkcm_indepTest()] together with
#' [gkcm_suffStat()] to provide this test in the interface required by packages such as
#' `pcalg`.
#'
#' @param x A vector of observations for the first variable. Must be an atomic vector
#'   (numeric/integer/logical/factor/character), contain no missing values, have at least
#'   10 observations, and must not be constant. Numeric values must be finite (no `Inf`,
#'   `-Inf`, or `NaN`).
#' @param y A vector of observations for the second variable. Same requirements as `x`,
#'   and must have the same length as `x`.
#' @param Z Conditioning variables. Either `NULL` (marginal test) or a `data.frame`,
#'   numeric/integer matrix, or atomic vector (treated as a single conditioning variable).
#'   Must contain no missing values, have at least 10 observations, and must not contain
#'   constant columns/values. Numeric values must be finite. If `Z` is not `NULL`, its number
#'   of observations must match `length(x)` and `length(y)`. If a 0-column `data.frame` or
#'   matrix is supplied, it is interpreted as `Z = NULL` (with a warning).
#' @param num_threads Either `NULL` or a single positive integer. The resulting
#'   number of threads is passed to `drf` and `ranger`. If `num_threads = NULL`,
#'   half of the available cores are used (minimum 1). Set `num_threads = 1`
#'   to disable parallelization.
#'
#' @return A list with components:
#' \describe{
#'   \item{t_stat}{Non-negative numeric scalar.}
#'   \item{p_val}{Numeric scalar in `[0, 1]`.}
#' }
#'
#' @family GKCM conditional independence tests
#' @seealso [gkcm_indepTest()], [gkcm_suffStat()]
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' z <- rnorm(n)
#' x <- z + rnorm(n, sd = 0.5)
#' y <- z + rnorm(n, sd = 0.5)
#'
#' gkcm(x, y, Z = z)            # conditional
#' gkcm(x, y, Z = NULL)         # marginal
#' gkcm(x, y, Z = data.frame()) # treated as marginal (warns)
#'
#' g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
#' gkcm(x, y, Z = data.frame(z = z, g = g))
#'
#' @export
gkcm <- function(x, y, Z, num_threads = NULL){

  Z <- .normalize_cond(Z)
  x <- .normalize_var(x, Z)
  y <- .normalize_var(y, Z)
  .assert_n(x, y)
  .assert_num_threads(num_threads)

  n_threads <- if (is.null(num_threads)) {
    # Use half of available cores (min. 1)
    max(1, floor(parallelly::availableCores() / 2))
  } else {
    num_threads
  }

  K_r <- .get_residual_gram(x, Z, num_threads = n_threads) |> .center_gram()
  L_r <- .get_residual_gram(y, Z, num_threads = n_threads) |> .center_gram()

  R <- K_r * L_r
  n <- ncol(R)
  sum_R <-

  # Compute test statistic
  t_stat <- sum(R) / n

  # Estimate eigenvalues
  H <- .center_gram(R)/(n - 1)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values

  p_val <- .compute_p_val(t_stat, ev)

  list(t_stat = t_stat, p_val = p_val)

}

#' GKCM conditional independence test for causal discovery
#'
#' Computes a p-value for testing whether variable `x` is conditionally independent
#' of variable `y` given the set of variables `S`. This function is intended to be
#' supplied as `indepTest` in constraint-based causal discovery algorithms.
#'
#' @details
#' This function is called as `indepTest(x, y, S, suffStat)`, where `x` and `y` are
#' integer indices of the tested variables and `S` is an integer vector of conditioning
#' indices (possibly empty). All indices refer to column positions in `suffStat$df`.
#' The kernel embedding is determined by the column types in `suffStat$df` (Gaussian,
#' \eqn{\sigma = 1}, for numeric variables; Dirac for categorical variables).
#'
#' `suffStat` must be created by [gkcm_suffStat()] and is passed through by algorithms
#' such as [pcalg::pc()], [pcalg::fci()], [tpc::tpc()], and [causalDisco::tfci()].
#' For performance, `suffStat` is not re-validated on each call.
#'
#' @param x,y Integer scalar indices of the variables to test (column positions in `suffStat$df`).
#' @param S Integer vector of conditioning set indices. May be empty (`integer(0)`), in which
#'   case a marginal (unconditional) independence test is performed.
#' @param suffStat A list of sufficient statistics created by [gkcm_suffStat()].
#' @param num_threads Either `NULL` or a single positive integer. The resulting
#'   number of threads is passed to `drf` and `ranger`. If `num_threads = NULL`,
#'   half of the available cores are used (minimum 1). Set `num_threads = 1`
#'   to disable parallelization.
#'
#' @return Numeric scalar in `[0, 1]`.
#'
#' @family GKCM conditional independence tests
#' @seealso [gkcm()], [gkcm_suffStat()], [pcalg::pc()], [pcalg::fci()], [tpc::tpc()], [causalDisco::tfci()]
#'
#' @examples
#' if (requireNamespace("pcalg", quietly = TRUE)) {
#'   df <- airquality[complete.cases(airquality), ]
#'   suffStat <- gkcm_suffStat(df)
#'   g <- pcalg::fci(
#'     suffStat = suffStat,
#'     indepTest = gkcm_indepTest,
#'     alpha = 0.2,
#'     labels = names(df)
#'   )
#' }
#'
#' @export
gkcm_indepTest <- function(x, y, S, suffStat, num_threads = NULL){

  n_threads <- if (is.null(num_threads)) {
    # Use half of available cores (min. 1)
    max(1, floor(parallelly::availableCores() / 2))
  } else {
    num_threads
  }

  if (length(S) == 0) {

    K_r <- .get_residual_gram(
      suffStat[[1]][, x], NULL,
      suffStat[[2]][[x]]
    ) |> .center_gram()

    L_r <- .get_residual_gram(
      suffStat[[1]][, y], NULL,
      suffStat[[2]][[y]]
    ) |> .center_gram()

  } else {

    K_r <- .get_residual_gram(
      suffStat[[1]][, x],
      suffStat[[1]][, S, drop = F],
      suffStat[[2]][[x]],
      num_threads = n_threads
    ) |> .center_gram()

    L_r <- .get_residual_gram(
      suffStat[[1]][, y],
      suffStat[[1]][, S, drop = F],
      suffStat[[2]][[y]],
      num_threads = n_threads
    ) |> .center_gram()

  }

  R <- K_r * L_r
  n <- ncol(R)

  # Compute test statistic
  t_stat <- sum(R) / n

  # Estimate eigenvalues
  H <- .center_gram(R)/(n - 1)
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values

  .compute_p_val(t_stat, ev)

}

