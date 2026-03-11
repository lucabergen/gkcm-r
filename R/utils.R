
# Checks the format of x/y and returns a centered numeric or factor vector
#' @noRd
.normalize_var <- function(x, Z = NULL) {

  n <- if (is.null(Z)) NULL else nrow(Z)
  .assert_var(x, n = n, .var.name = deparse(substitute(x)))

  if (is.numeric(x)) {
    return(x - mean(x))
  }
  factor(x)

}

# Checks the format of Z and returns either NULL or coerces to a df with numeric
# and factor columns
#' @noRd
.normalize_cond <- function(Z) {

  # Treat a 0-column conditioning object as "no conditioning"
  if (!is.null(Z) && (is.data.frame(Z) || is.matrix(Z)) && ncol(Z) == 0L) {
    warning(
      "`Z` has 0 columns and is interpreted as `Z = NULL`; performing a marginal independence test.",
      call. = FALSE
    )
    return(NULL)
  }

  .assert_cond(Z)

  if (is.null(Z)) return(NULL)

  # Coerce to data.frame
  df <- if (is.data.frame(Z)) {
    Z
  } else if (is.matrix(Z)) {
    as.data.frame(Z, stringsAsFactors = FALSE)
  } else {
    # vector -> 1-column df
    data.frame(Z = Z, check.names = FALSE, stringsAsFactors = FALSE)
  }

  # Transform character + logical columns to factors
  idx <- vapply(df, function(col) is.character(col) || is.logical(col), logical(1))
  df[idx] <- lapply(df[idx], factor)

  df
}

# Checks the format of df and returns a df with numeric and factor columns
#' @noRd
.normalize_df <- function(df) {

  .assert_df(df)

  # Transform character + logical columns to factors
  idx <- vapply(df, function(col) is.character(col) || is.logical(col), logical(1))
  df[idx] <- lapply(df[idx], factor)

  df

}

#' @noRd
.gaussian_gram <- function(x, sigma = 1) {

  x <- as.vector(x)
  exp(-outer(x, x, "-")^2 / (2 * sigma^2))

}

#' @noRd
.dirac_gram <- function(x) {

  outer(as.integer(x), as.integer(x), FUN = "==") * 1.0

}

#' @noRd
.get_gram <- function(x) {

  if (is.factor(x)) {
    .dirac_gram(x)
  } else {
    .gaussian_gram(x)
  }

}

#' @noRd
.center_gram <- function(K, symmetrize = TRUE) {

  r <- rowMeans(K)
  m <- mean(K)

  # Center: Kc = K - r*1' - 1*r' + m
  Kc <- sweep(K, 1, r, FUN = "-")
  Kc <- sweep(Kc, 2, r, FUN = "-")
  Kc <- Kc + m

  # Symmetrize
  (Kc + t(Kc)) / 2

}

#' @noRd
.one_hot_encode <- function(x){

  classes <- levels(x)
  n <- length(x)
  k <- length(classes)

  OH <- matrix(0L, nrow = n, ncol = k)
  colnames(OH) <- classes

  OH[cbind(seq_len(n), as.integer(x))] <- 1L

  OH

}

#' @noRd
.get_feat_repr <- function(x) {

  if (is.factor(x)) {
    .one_hot_encode(x)
  } else {
    .gaussian_gram(x)
  }

}

#' GKCM sufficient statistics for causal discovery
#'
#' Prepares sufficient statistics for [gkcm_indepTest()] in constraint-based causal
#' discovery algorithms (e.g., in `pcalg`).
#'
#' @details
#' Returns a `suffStat` object containing a normalized copy of `df` and per-variable
#' feature representations, using a Gaussian kernel (bandwidth \eqn{\sigma = 1}) for
#' numeric variables and a Dirac kernel (via one-hot encoding) for categorical variables.
#'
#' @param df A data frame with one column per variable and one row per observation.
#'   Must have at least 10 rows and at least 2 columns, contain no missing values,
#'   and only use column types numeric/integer/logical/factor/character. Numeric
#'   values must be finite (no `Inf`, `-Inf`, `NaN`), and columns must not be constant.
#'
#' @return A list with components:
#' \describe{
#'   \item{df}{The normalized input data. Variable indices used by causal discovery
#'     algorithms refer to the column positions of this data frame.}
#'   \item{feat_repr_list}{A list of per-variable feature representations used by
#'     [gkcm_indepTest()].}
#' }
#'
#' @family GKCM conditional independence tests
#' @seealso [gkcm()], [gkcm_indepTest()], [pcalg::pc()], [pcalg::fci()]
#'
#' @examples
#' df <- data.frame(x = rnorm(100), y = rnorm(100), z = rnorm(100))
#' suffStat <- gkcm_suffStat(df)
#' str(suffStat)
#'
#' df <- airquality[complete.cases(airquality), ]
#' suffStat <- gkcm_suffStat(df)
#' str(suffStat)
#'
#' @export
gkcm_suffStat <- function(df){

  .assert_df(df)

  df <- .normalize_df(df)
  feat_repr_list <- lapply(df, .get_feat_repr)

  list(df = df, feat_repr_list = feat_repr_list)

}

#' @noRd
.compute_drf_residuals <- function(K, W) {

  n <- nrow(K)

  # Compute (I - W)K(I - W)' =
  # K - WK - KW' + WKW' = K - WK - (WK)' + WKW' (K symmetric)

  K <- Matrix::Matrix(K, sparse = FALSE)

  WK   <- W %*% K
  WKWt <- WK %*% Matrix::t(W)
  R <- K - WK - Matrix::t(WK) + WKWt
  R <- Matrix::forceSymmetric(R, uplo = "U")

  as.matrix(R)

  # WK <- W %*% K
  # R  <- K - WK - t(WK) + WK %*% t(W)
  #
  # # Symmetrize
  # (R + t(R)) / 2

}

#' @noRd
.compute_p_val <- function(t_stat, ev) {

  # Filter nonzero eigenvalues
  tol <- 100 * .Machine$double.eps * max(1, max(abs(ev)))
  ev[ev < 0 & abs(ev) <= tol] <- 0
  ev <- ev[ev > tol]

  if (length(ev) == 0) {
    p_val <- if (t_stat <= tol) {1} else {0}
  } else if (t_stat <= 0) {
    p_val <- 1
  } else {
    dvs <- CompQuadForm::davies(t_stat, lambda = ev)
    if (dvs$ifault == 0 && is.finite(dvs$Qq) && dvs$Qq >= 0 && dvs$Qq <= 1) {
      p_val <- dvs$Qq
    } else {
      imh <- CompQuadForm::imhof(t_stat, lambda = ev)
      p_val <- min(1, max(0, imh$Qq))
    }
  }

  p_val
}
