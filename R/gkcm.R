#' GKCM conditional independence test
#'
#' Tests whether the random variables X and Y are conditionally independent
#' given Z using the Generalised Kernel Covariance Measure (GKCM). Their
#' observed data are supplied through arguments `x` and `y`, and the
#' conditioning data through `z`. Supports numeric and categorical variables,
#' custom feature or Gram-matrix representations, and RKHS-valued
#' regression with `okranger`.
#'
#' @details
#' By default, numeric vectors supplied through `x` and `y` are transformed
#' using random Fourier features (RFFs, `num_RFFs = 20`). These approximate the
#' feature map of the Gaussian kernel with unit lengthscale without constructing
#' Gram matrices. Set `num_RFFs = NULL` to use the Gaussian kernel instead.
#' Categorical variables supplied through `x` or `y` must be encoded as factors.
#' They use the Dirac kernel, whose feature map can be represented by a one-hot
#' encoding of the categories. Missing values are not allowed. Exclude or impute
#' these before testing.
#'
#' In conditional tests, the feature-mapped variables are regressed on the
#' conditioning variables using random forests. `gkcm()` uses `okranger` for
#' feature and Gram matrices and `ranger` for factor representations. The
#' default values for tuning parameters are `mtry = ncol(z)`, `num.trees = 200`,
#' `min.node.size = 5`, `max.depth = 0` (no depth limit), and
#' `approx.tolerance = 0.001` for feature and Gram matrices and
#' `mtry = ncol(z)`, `num.trees = 200`, `min.node.size = 10`, and
#' `max.depth = 0` for factor representations. These defaults can be overridden
#' using `ranger_args_x` or `ranger_args_y`; unspecified defaults remain in
#' effect. `mtry` may be a function. For each conditional test, the regression
#' backend calls it with `ncol(z)` as its only argument. It must return a single
#' whole number that is at least 1 and at most `ncol(z)`. This allows `mtry` to
#' adapt to the number of conditioning variables. `num_threads` sets the number
#' of threads used by the random forests. With `z = NULL`, GKCM performs a
#' marginal independence test and skips regression.
#'
#' Custom feature or Gram matrices can be supplied for X and Y (`Phi` or `K` and
#' `Psi` or `L`, respectively). Gram matrices must be symmetric and positive
#' semidefinite; positive semidefiniteness is not checked. Since feature and
#' Gram matrices use a different regression backend and tuning defaults than
#' factor representations, we recommend supplying custom matrices only for
#' numeric variables. For categorical variables, supply factor vectors through
#' `x` and `y`. If multiple representations are given for a variable, the
#' function warns and uses the feature matrix in preference to the Gram matrix,
#' and the Gram matrix in preference to the data. Supply feature matrices
#' instead of Gram matrices when the features are not high-dimensional.
#'
#' For causal discovery workflows (e.g., PC/FCI), use [gkcm_indepTest()]
#' together with [gkcm_suffStat()].
#'
#' @param x Numeric or integer vector for numeric X, or a factor vector if X is
#'   categorical. Alternatively, supply `Phi` or `K`.
#' @param y Numeric or integer vector for numeric Y, or a factor vector if Y is
#'   categorical. Alternatively, supply `Psi` or `L`.
#' @param z `NULL` for a marginal test, or a data frame of conditioning
#'   variables with at least one column. Columns must be numeric or integer for
#'   continuous variables or factors for categorical variables. Missing values
#'   are not allowed.
#' @param num_RFFs Even positive integer giving the number of random Fourier
#'   features used for numeric vector inputs `x` and `y`. Defaults to `20`.
#'   Set to `NULL` to use the Gaussian kernel directly.
#' @param num_threads Number of threads used by the regression forests. If
#'   `NULL`, half of the available CPU cores are used (at least one). Set to `1`
#'   to use a single thread.
#' @param seed Optional integer seed for random Fourier features and regression
#'   fits. If `NULL`, random seeds are drawn.
#' @param K Optional Gram-matrix representation for X. One row and
#'   column per observation.
#' @param L Optional Gram-matrix representation for Y. One row and
#'   column per observation.
#' @param Phi Optional feature-matrix representation for X. One row per
#'   observation and at least one feature column.
#' @param Psi Optional feature-matrix representation for Y. One row per
#'   observation and at least one feature column.
#' @param ranger_args_x If not `NULL`, a non-empty named list of supported
#'   tuning arguments for regressing the features of X on Z.
#' @param ranger_args_y If not `NULL`, a non-empty named list of supported
#'   tuning arguments for regressing the features of Y on Z.
#'
#' @return A list with components:
#' \describe{
#'   \item{statistic}{Non-negative numeric scalar.}
#'   \item{p_value}{Numeric scalar in `[0, 1]`.}
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
#' gkcm(x, y, z = data.frame(z = z)) # conditional
#' gkcm(x, y, z = NULL)              # marginal
#'
#' g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
#' gkcm(x, y, z = data.frame(z = z, g = g))
#'
#' # Pass categorical tested variables as factors:
#' gkcm(x, g, z = data.frame(z = z))
#'
#' # User-supplied representations are also supported:
#' Phi <- cbind(cos(x), sin(x))
#' gkcm(Phi = Phi, y = y, z = data.frame(z = z))
#'
#' @export
gkcm <- function(x = NULL, y = NULL, z = NULL,
                 num_RFFs = 20, num_threads = NULL, seed = NULL,
                 K = NULL, L = NULL, Phi = NULL, Psi = NULL,
                 ranger_args_x = NULL, ranger_args_y = NULL){

  seed <- .validate_seed(seed)
  num_RFFs <- .validate_num_RFFs(num_RFFs)
  num_threads <- .validate_num_threads(num_threads)

  x_rep <- .choose_representation(
    args = list(x = x, K = K, Phi = Phi),
    var = "x"
  )
  x_rep <- .validate_representation(x_rep, "x")

  y_rep <- .choose_representation(
    args = list(y = y, L = L, Psi = Psi),
    var = "y"
  )
  y_rep <- .validate_representation(y_rep, "y")

  n <- NROW(x_rep$value)

  if (NROW(y_rep$value) != n) {
    stop("The data must have the same number of observations.",
         call. = FALSE)
  }

  z <- .validate_z(z, n)

  if (.representation_is_constant(x_rep)) {
    warning("The representation of x is constant; returning statistic = 0 and p_value = 1.")
    return(list(statistic = 0, p_value = 1))
  }

  if (.representation_is_constant(y_rep)) {
    warning("The representation of y is constant; returning statistic = 0 and p_value = 1.")
    return(list(statistic = 0, p_value = 1))
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  seed_x <- sample.int(.Machine$integer.max, 1)
  seed_y <- sample.int(.Machine$integer.max, 1)

  # If numeric data was given, save Gaussian Gram matrix or Random Fourier Features
  if (x_rep$type == "data") {
    x_rep <- .transform_data(x_rep, num_RFFs)
  }

  if (y_rep$type == "data") {
    y_rep <- .transform_data(y_rep, num_RFFs)
  }

  ranger_args_x <- .validate_ranger_args(ranger_args_x, "x", x_rep$type)
  ranger_args_y <- .validate_ranger_args(ranger_args_y, "y", y_rep$type)

  # Check feature dimensions for explicit representations
  if (x_rep$type != "gram" && y_rep$type != "gram") {

    if (x_rep$type == "features") {
      x_dim <- ncol(x_rep$value)
    } else {
      x_dim <- nlevels(x_rep$value)
    }

    if (y_rep$type == "features") {
      y_dim <- ncol(y_rep$value)
    } else {
      y_dim <- nlevels(y_rep$value)
    }

    if (x_dim * y_dim > 20^2 && x_dim * y_dim >= n) {
      warning("The feature dimensions are large, so GKCM may be slow. Consider using Gram matrices as inputs.")
    }
  }

  if (is.null(z)) {

    eps_rep <- .center_representation(x_rep)
    xi_rep  <- .center_representation(y_rep)

  } else {

    eps_rep <- .regress_representation(
      x_rep, z,
      num_threads = num_threads,
      seed = seed_x,
      ranger_args = ranger_args_x
    )
    xi_rep  <- .regress_representation(
      y_rep, z,
      num_threads = num_threads,
      seed = seed_y,
      ranger_args = ranger_args_y
    )

    eps_rep <- .center_representation(eps_rep)
    xi_rep  <- .center_representation(xi_rep)

  }

  components <- .compute_test_components(eps_rep, xi_rep)

  p_value <- .compute_p_value(components)

  return(list(statistic = components$statistic, p_value = p_value))

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
#' When values in `suffStat$df` are missing, `indepTest` uses testwise deletion,
#' where rows with a missing value in `x`, `y`, or any variable in `S` are
#' removed without a warning.
#'
#' `suffStat` must be created by [gkcm_suffStat()] and is passed through by algorithms
#' such as [pcalg::pc()], [pcalg::fci()], [tpc::tpc()], and [causalDisco::tfci()].
#' For performance, `suffStat` is not validated by [gkcm_suffStat()].
#'
#' @param x,y Integer scalar indices of the variables to test (column positions in `suffStat$df`).
#' @param S Integer vector of conditioning set indices. May be empty (`integer(0)`), in which
#'   case a marginal independence test is performed.
#' @param suffStat A list of sufficient statistics created by [gkcm_suffStat()].
#'
#' @return p-value in `[0, 1]`.
#'
#' @family GKCM conditional independence tests
#' @seealso [gkcm()], [gkcm_suffStat()], [pcalg::pc()], [pcalg::fci()], [tpc::tpc()], [causalDisco::tfci()]
#'
#' @examples
#' if (requireNamespace("pcalg", quietly = TRUE)) {
#'   df <- iris
#'   suffStat <- gkcm_suffStat(df)
#'   g <- pcalg::fci(
#'     suffStat = suffStat,
#'     indepTest = gkcm_indepTest,
#'     alpha = 0.1,
#'     labels = names(df)
#'   )
#' }
#'
#' @export
gkcm_indepTest <- function(x, y, S, suffStat) {

  if (!inherits(suffStat, "gkcm_suffStat")) {
    stop(
      "'suffStat' must be created by gkcm_suffStat().",
      call. = FALSE
    )
  }

  # Accept NULL as empty conditioning set and sort predictors
  if (is.null(S)) {
    S <- integer(0)
  } else {
    S <- sort(S)
  }

  # Testwise deletion: select only complete cases for x, y, and S
  used <- c(x, y, S)

  rows <- which(stats::complete.cases(suffStat$df[, used, drop = FALSE]))

  n <- length(rows)

  if (n < 2) {
    stop(
      "At least 2 complete observations are required for this test.",
      call. = FALSE
    )
  }

  x_rep <- suffStat$representations[[names(suffStat$df)[x]]]
  y_rep <- suffStat$representations[[names(suffStat$df)[y]]]

  # Subset the representations using the same observations
  if (x_rep$type == "gram") {
    x_rep$value <- x_rep$value[rows, rows, drop = FALSE]
  } else if (x_rep$type == "features") {
    x_rep$value <- x_rep$value[rows, , drop = FALSE]
  } else {
    x_rep$value <- droplevels(x_rep$value[rows])
  }

  if (y_rep$type == "gram") {
    y_rep$value <- y_rep$value[rows, rows, drop = FALSE]
  } else if (y_rep$type == "features") {
    y_rep$value <- y_rep$value[rows, , drop = FALSE]
  } else {
    y_rep$value <- droplevels(y_rep$value[rows])
  }

  if (length(S) == 0) {
    z <- NULL
  } else {
    z <- suffStat$df[rows, S, drop = FALSE]
  }

  z <- .validate_z(z, n)

  ranger_args <- suffStat$ranger_args
  if (is.null(ranger_args)) {
    ranger_args <- list()
  }
  ranger_args_x <- ranger_args[[names(suffStat$df)[x]]]
  ranger_args_y <- ranger_args[[names(suffStat$df)[y]]]

  # Check for constant representations after testwise deletion
  if (.representation_is_constant(x_rep)) {
    warning(
      sprintf(
        "The representation of x is constant for the test x = %d, y = %d, S = {%s}; returning p-value 1.",
        x, y, paste(S, collapse = ", ")
      ),
      call. = FALSE
    )
    return(1)
  }

  if (.representation_is_constant(y_rep)) {
    warning(
      sprintf(
        "The representation of y is constant for the test x = %d, y = %d, S = {%s}; returning p-value 1.",
        x, y, paste(S, collapse = ", ")
      ),
      call. = FALSE
    )
    return(1)
  }

  if (is.null(z)) {

    eps_rep <- .center_representation(x_rep)
    xi_rep  <- .center_representation(y_rep)

  } else {

    if (is.null(suffStat$seed)) {
      seed_x <- sample.int(.Machine$integer.max, 1)
      seed_y <- sample.int(.Machine$integer.max, 1)
    } else {
      # Derive deterministic seeds
      seed_x <- .derive_regression_seed(suffStat$seed, x, y, S, target = x)
      seed_y <- .derive_regression_seed(suffStat$seed, x, y, S, target = y)
    }

    eps_rep <- .regress_representation(
      x_rep, z,
      num_threads = suffStat$num_threads,
      seed = seed_x,
      ranger_args = ranger_args_x
    )
    xi_rep  <- .regress_representation(
      y_rep, z,
      num_threads = suffStat$num_threads,
      seed = seed_y,
      ranger_args = ranger_args_y
    )

    eps_rep <- .center_representation(eps_rep)
    xi_rep  <- .center_representation(xi_rep)

  }

  components <- .compute_test_components(eps_rep, xi_rep)

  .compute_p_value(components)
}

