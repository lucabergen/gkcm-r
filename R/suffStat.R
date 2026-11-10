
#' Prepare GKCM independence tests for causal discovery
#'
#' Creates a `gkcm_suffStat` object for use with [gkcm_indepTest()] in
#' constraint-based causal discovery. It stores the data, representations
#' of the variables, and settings used by the independence tests.
#'
#' @details
#' By default, numeric columns are transformed using random Fourier features
#' (RFFs, `num_RFFs = 20`). These approximate the feature map of the Gaussian
#' kernel with unit lengthscale without constructing Gram matrices. Set
#' `num_RFFs = NULL` to use the Gaussian kernel instead. Categorical variables
#' must be encoded as factors. They use the Dirac kernel, whose feature map can
#' be represented by a one-hot encoding of the categories.
#'
#' Custom feature or Gram matrices can be supplied for variables in `df`. These
#' take precedence over automatically generated representations, independently
#' of `num_RFFs`. Gram matrices must be symmetric and positive semidefinite;
#' positive semidefiniteness is not checked. Since feature and Gram matrices use
#' a different regression backend and tuning parameter defaults than factor
#' representations (see below), we recommend supplying custom matrices only for
#' numeric variables. This preserves the backend and parameter defaults for
#' factor representations of categorical variables.
#'
#' In conditional tests, the feature-mapped variables are regressed on the
#'  variables using random forests. [gkcm_indepTest()] uses
#' `okranger` for feature and Gram matrices and `ranger` for factor
#' representations. The default values for tuning parameters are
#' `mtry = length(S)`, `num.trees = 200`, `min.node.size = 5`, `max.depth = 0`
#' (no depth limit), and `approx.tolerance = 0.001` for feature and Gram
#' matrices and `mtry = length(S)`, `num.trees = 200`, `min.node.size = 10`, and
#' `max.depth = 0` for factor representations, where `S` is the
#' conditioning set. These defaults can be overridden per column using
#' `ranger_args`; unspecified defaults remain in effect. `mtry` may be a
#' function. For each conditional test, the regression backend calls it with
#' `length(S)` as its only argument. It must return a single integer that is at
#' least 1 and at most `length(S)`. This allows `mtry` to adapt to the size of
#' the conditioning set. `num_threads` sets the number of threads used by the
#' random forests. Marginal tests skip regression.
#'
#' Missing values are allowed in `df`. [gkcm_indepTest()] applies testwise
#' deletion without a warning for each test, using only observations complete
#' for the tested and conditioning variables. In user-supplied representations,
#' missing values must match those in the corresponding column of `df`, i.e.,
#' `NA` rows in feature matrices, and `NA` rows and columns in Gram matrices.
#'
#' @param df A data frame with at least two rows and unique, non-empty column
#'   names. Columns must be numeric/integer or factors for categorical variables.
#'   Missing values are allowed, but no column may be entirely missing; observed
#'   numeric values must be finite.
#' @param representations Optional named list of representations, indexed by
#'   column names in `df`. The list may be incomplete or empty; omitted
#'   representations are generated from `df`. Each entry is a list with `type`
#'   and `value`: `"features"` for a numeric feature matrix, `"gram"` for a
#'   numeric Gram matrix, or `"data"` for a factor vector. Feature matrices have
#'   one row per observation; Gram matrices have one row and column per observation.
#'   See Details.
#' @param num_RFFs Even positive integer giving the number of random Fourier
#'   features used for automatically generated numeric representations; defaults
#'   to `20`. Set to `NULL` to use the Gaussian kernel directly.
#' @param num_threads Number of threads used by the regression forests. If `NULL`,
#'   half of the available CPU cores are used (at least one). Set to `1` to use
#'   a single thread.
#' @param seed Optional integer seed for generating random Fourier features and
#'   regression fits. If `NULL`, a random seed is drawn.
#' @param ranger_args Optional named list of regression argument lists, indexed
#'   by the column names of `df`. The list may be incomplete or empty. Each
#'   entry must be a non-empty named list of supported tuning arguments.
#'   Omit a column to use its defaults. See Details.
#'
#' @return An object of class `"gkcm_suffStat"` containing the data frame,
#'   representations, per-column regression arguments, and test settings.
#'
#' @family GKCM conditional independence tests
#' @seealso [gkcm()], [gkcm_indepTest()]
#'
#' @examples
#' df <- iris
#' suffStat <- gkcm_suffStat(df, num_RFFs = 10, seed = 1)
#'
#' # Test Sepal.Length and Sepal.Width given Petal.Length:
#' gkcm_indepTest(1, 2, 3, suffStat)
#'
#' # Adapt mtry to the conditioning set for a numeric and a factor variable:
#' mtry_rule <- function(p) min(2, p)
#' suffStat <- gkcm_suffStat(
#'   df, num_RFFs = 10, seed = 1,
#'   ranger_args = list(
#'     Sepal.Length = list(mtry = mtry_rule),
#'     Species = list(mtry = mtry_rule)
#'   )
#' )
#' gkcm_indepTest(1, 5, c(2, 3), suffStat)
#'
#' # Supply a feature representation for a numeric variable:
#' Phi <- cbind(cos(df$Sepal.Length), sin(df$Sepal.Length))
#' suffStat <- gkcm_suffStat(
#'   df, num_RFFs = 10,
#'   representations = list(Sepal.Length = list(type = "features", value = Phi))
#' )
#'
#' @export
gkcm_suffStat <- function(df, representations = NULL, num_RFFs = 20,
                          num_threads = NULL, seed = NULL, ranger_args = NULL) {

  .validate_suffStat_df(df)

  num_RFFs <- .validate_num_RFFs(num_RFFs)
  num_threads <- .validate_num_threads(num_threads)
  seed <- .validate_seed(seed)

  if (is.null(representations)) {
    representations <- list()
  }

  if (!is.list(representations)) {
    stop(
      "'representations' must be NULL or a named list.",
      call. = FALSE
    )
  }

  column_names <- names(df)
  supplied_names <- names(representations)

  if (length(representations) > 0L) {

    unknown_names <- setdiff(supplied_names, column_names)

    if (length(unknown_names) > 0) {
      stop(
        sprintf(
          "Unknown column names in 'representations': %s.",
          paste(unknown_names, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    # Validate supplied entries before generating missing representations
    .validate_suffStat_representations(
      representations = representations,
      df = df[, supplied_names, drop = FALSE]
    )
  } else {
    supplied_names <- character(0)
  }

  missing_names <- setdiff(column_names, supplied_names)

  if (!is.null(num_RFFs) && num_RFFs >= max(21, sqrt(nrow(df)))) {
    warning("'num_RFFs' is large relative to the sample size. Set 'num_RFFs' < sqrt(n) to speed up causal discovery.")
  }

  if (length(missing_names) > 0) {
    missing_df <- df[, missing_names, drop = FALSE]

    generated <- .generate_suffStat_representations(
      df = missing_df,
      num_RFFs = num_RFFs,
      seed = seed
    )

    .validate_suffStat_representations(
      representations = generated,
      df = missing_df
    )

    representations <- c(representations, generated)
  }

  # Store the complete list in data-frame column order
  representations <- representations[column_names]

  ranger_args <- .validate_suffStat_ranger_args(
    ranger_args = ranger_args,
    df = df,
    representations = representations
  )

  dimensions <- sort(vapply(
    representations,
    function(rep) {
      if (rep$type == "gram") {
        0
      } else if (rep$type == "features") {
        ncol(rep$value)
      } else {
        nlevels(droplevels(rep$value[!is.na(rep$value)]))
      }
    },
    numeric(1)
  ), decreasing = TRUE)

  if (length(dimensions) >= 2 &&
      dimensions[1] * dimensions[2] > 20^2 &&
      dimensions[1] * dimensions[2] >= nrow(df)) {
    warning(
      "Some feature dimensions are large relative to the sample size; GKCM tests may be slow.",
      call. = FALSE
    )
  }

  if (any(vapply(df, anyNA, logical(1)))) {
    warning(
      paste0(
        "Missing values are present in `df`. `gkcm_indepTest()` will apply ",
        "testwise deletion without a warning for each test."
      ),
      call. = FALSE
    )
  }

  .new_gkcm_suffStat(
    df = df,
    representations = representations,
    ranger_args = ranger_args,
    num_RFFs = num_RFFs,
    num_threads = num_threads,
    seed = seed
  )

}

#' Constructor for gkcm_suffStat objects
#'
#' @noRd
.new_gkcm_suffStat <- function(df, representations, ranger_args, num_RFFs,
                               num_threads, seed) {

  stopifnot(
    is.data.frame(df),
    is.list(representations)
  )

  structure(
    list(
      df = df,
      representations = representations,
      ranger_args = ranger_args,
      num_RFFs = num_RFFs,
      num_threads = num_threads,
      seed = seed
    ),
    class = "gkcm_suffStat"
  )
}

#' Generate representations for GKCM sufficient statistics; assumes that
#' df, num_RFFs and seed have already been validated.
#'
#' @noRd
.generate_suffStat_representations <- function(df, num_RFFs, seed) {

  n <- nrow(df)

  representations <- vector("list", ncol(df))
  names(representations) <- names(df)

  has_numeric <- any(vapply(
    df,
    function(column) {
      is.numeric(column) && !is.factor(column)
    },
    logical(1)
  ))

  # Set the seed for generating numeric representations
  if (has_numeric && !is.null(seed)) {
    set.seed(seed)
  }

  for (j in seq_along(df)) {
    column <- df[[j]]
    observed <- !is.na(column)

    if (is.numeric(column) && !is.factor(column)) {
      observed_values <- column[observed]

      if (is.null(num_RFFs)) {
        observed_rep <- .gaussian_kernel(observed_values)

        value <- matrix(
          NA_real_,
          nrow = n,
          ncol = n
        )

        value[observed, observed] <- observed_rep

        representations[[j]] <- list(
          type = "gram",
          value = value
        )

      } else {
        observed_rep <- .rff(observed_values, num_RFFs)

        value <- matrix(
          NA_real_,
          nrow = n,
          ncol = ncol(observed_rep)
        )

        value[observed, ] <- observed_rep

        representations[[j]] <- list(
          type = "features",
          value = value
        )
      }

    } else {
      if (is.factor(column)) {
        # Unused factor levels are dropped in gkcm_indepTest
        value <- column
      } else {
        value <- factor(column)
      }

      representations[[j]] <- list(
        type = "data",
        value = value
      )
    }
  }

  representations
}
