
#' @noRd
.validate_representation <- function(rep, var) {
  value <- rep$value

  if (rep$type == "data") {
    checkmate::assert_true(
      is.atomic(value) && is.null(dim(value)) && !is.data.frame(value),
      .var.name = var
    )

    if (is.numeric(value) && !is.complex(value)) {
      checkmate::assert_numeric(
        value,
        any.missing = FALSE,
        finite = TRUE,
        min.len = 2,
        .var.name = var
      )
    } else if (is.factor(value)) {
      checkmate::assert_factor(
        value,
        any.missing = FALSE,
        min.len = 2,
        .var.name = var
      )
      value <- droplevels(value)
    } else if (is.character(value)) {
      checkmate::assert_character(
        value,
        any.missing = FALSE,
        min.len = 2,
        .var.name = var
      )
      value <- factor(value)
    } else if (is.logical(value)) {
      checkmate::assert_logical(
        value,
        any.missing = FALSE,
        min.len = 2,
        .var.name = var
      )
      value <- factor(value)
    } else {
      stop(
        sprintf(
          "'%s' must be a numeric, factor, character, or logical vector.",
          var
        ),
        call. = FALSE
      )
    }

  } else if (rep$type == "features") {
    checkmate::assert_matrix(
      value,
      mode = "numeric",
      any.missing = FALSE,
      min.rows = 2,
      min.cols = 1,
      .var.name = var
    )
    checkmate::assert_true(!is.complex(value), .var.name = var)
    checkmate::assert_numeric(
      as.vector(value),
      any.missing = FALSE,
      finite = TRUE,
      .var.name = var
    )

  } else if (rep$type == "gram") {
    checkmate::assert_matrix(
      value,
      mode = "numeric",
      any.missing = FALSE,
      min.rows = 2,
      min.cols = 2,
      .var.name = var
    )
    checkmate::assert_true(!is.complex(value), .var.name = var)
    checkmate::assert_numeric(
      as.vector(value),
      any.missing = FALSE,
      finite = TRUE,
      .var.name = var
    )
    checkmate::assert_true(nrow(value) == ncol(value), .var.name = var)

    symmetric <- isTRUE(all.equal(
      value,
      t(value),
      tolerance = sqrt(.Machine$double.eps),
      check.attributes = FALSE
    ))

    if (!symmetric) {
      stop(sprintf("The Gram matrix for '%s' must be symmetric.", var),
           call. = FALSE)
    }

  } else {
    stop(sprintf("Unknown representation type for '%s'.", var),
         call. = FALSE)
  }

  rep$value <- value
  rep
}

.validate_z <- function(z, n) {

  if (is.null(z)) {
    return(NULL)
  }

  checkmate::assert_data_frame(z, min.cols = 1)
  checkmate::assert_true(nrow(z) == n, .var.name = "z")

  for (j in seq_along(z)) {
    column <- z[[j]]
    var_name <- sprintf("z[[%d]]", j)

    checkmate::assert_true(
      is.atomic(column) && is.null(dim(column)),
      .var.name = var_name
    )

    if (is.numeric(column) && !is.complex(column)) {
      checkmate::assert_numeric(
        column,
        any.missing = FALSE,
        finite = TRUE,
        .var.name = var_name
      )
    } else if (is.factor(column)) {
      checkmate::assert_factor(
        column,
        any.missing = FALSE,
        .var.name = var_name
      )
      z[[j]] <- droplevels(column)
    } else if (is.character(column)) {
      checkmate::assert_character(
        column,
        any.missing = FALSE,
        .var.name = var_name
      )
      z[[j]] <- factor(column)
    } else if (is.logical(column)) {
      checkmate::assert_logical(
        column,
        any.missing = FALSE,
        .var.name = var_name
      )
      z[[j]] <- factor(column)
    } else {
      stop(
        sprintf(
          "Column %d of 'z' must be numeric, factor, character, or logical.",
          j
        ),
        call. = FALSE
      )
    }
  }

  z
}

#' @noRd
.validate_seed <- function(seed) {

  if (is.null(seed)) {
    return(NULL)
  }

  checkmate::assert_number(
    seed,
    lower = -.Machine$integer.max,
    upper = .Machine$integer.max,
    finite = TRUE,
    .var.name = "seed"
  )
  checkmate::assert_true(seed == floor(seed), .var.name = "seed")

  as.integer(seed)
}

#' @noRd
.validate_num_RFFs <- function(num_RFFs) {

  if (is.null(num_RFFs)) {
    return(NULL)
  }

  checkmate::assert_number(
    num_RFFs,
    lower = 2,
    upper = .Machine$integer.max,
    finite = TRUE,
    .var.name = "num_RFFs"
  )
  checkmate::assert_true(
    num_RFFs == floor(num_RFFs),
    .var.name = "num_RFFs"
  )
  checkmate::assert_true(
    num_RFFs %% 2 == 0,
    .var.name = "num_RFFs"
  )

  as.integer(num_RFFs)
}

#' @noRd
.validate_num_threads <- function(num_threads) {

  if (is.null(num_threads)) {
    # Use half of available cores (min. 1)
    num_threads <- max(1, floor(parallelly::availableCores() / 2))
  }

  checkmate::assert_number(
    num_threads,
    lower = 1,
    upper = .Machine$integer.max,
    finite = TRUE,
    .var.name = "num_threads"
  )
  checkmate::assert_true(
    num_threads == floor(num_threads),
    .var.name = "num_threads"
  )

  as.integer(num_threads)
}

#' @noRd
.validate_ranger_args <- function(args, var, feature_type) {

  if (is.null(args)) {
    return(NULL)
  }

  checkmate::assert_list(
    args,
    min.len = 1L,
    names = "named",
    .var.name = var
  )
  checkmate::assert_true(!is.data.frame(args), .var.name = var)

  if (anyDuplicated(names(args))) {
    stop(sprintf("`ranger_args` for `%s` cannot contain duplicate names.", var),
         call. = FALSE)
  }

  allowed_args <- switch(
    feature_type,
    data = c("mtry", "num.trees", "min.node.size", "max.depth"),
    features = c("mtry", "num.trees", "min.node.size", "max.depth", "approx.tolerance"),
    gram = c("mtry", "num.trees", "min.node.size", "max.depth", "approx.tolerance")
  )
  unsupported_args <- setdiff(names(args), allowed_args)

  if (length(unsupported_args) > 0L) {
    stop(
      sprintf(
        paste0(
          "`ranger_args` for `%s` only supports: %s. Unsupported argument(s): %s."
        ),
        var,
        paste(allowed_args, collapse = ", "),
        paste(unsupported_args, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  args
}

#' Validate the data frame for gkcm_suffStat
#'
#' @noRd
.validate_suffStat_df <- function(df) {

  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame.", call. = FALSE)
  }

  if (nrow(df) < 2) {
    stop(
      "'df' must contain at least 2 observations.",
      call. = FALSE
    )
  }

  if (ncol(df) < 1) {
    stop(
      "'df' must contain at least one column.",
      call. = FALSE
    )
  }

  column_names <- names(df)

  if (is.null(column_names) ||
      anyNA(column_names) ||
      any(!nzchar(column_names)) ||
      anyDuplicated(column_names)) {
    stop(
      "'df' must have unique, non-empty column names.",
      call. = FALSE
    )
  }

  for (j in seq_along(df)) {
    column <- df[[j]]
    column_name <- column_names[[j]]

    if (!is.atomic(column) || !is.null(dim(column))) {
      stop(
        sprintf(
          "Column '%s' must be an atomic vector.",
          column_name
        ),
        call. = FALSE
      )
    }

    if (is.factor(column)) {
      # Factors are accepted as categorical variables.

    } else if (is.numeric(column) && !is.complex(column)) {
      if (any(is.nan(column)) || any(is.infinite(column))) {
        stop(
          sprintf(
            "Column '%s' must not contain NaN or infinite values.",
            column_name
          ),
          call. = FALSE
        )
      }

    } else if (!is.character(column) && !is.logical(column)) {
      stop(
        sprintf(
          "Column '%s' must be numeric, factor, character, or logical.",
          column_name
        ),
        call. = FALSE
      )
    }

    if (all(is.na(column))) {
      stop(
        sprintf(
          "Column '%s' must not be entirely missing.",
          column_name
        ),
        call. = FALSE
      )
    }
  }

  invisible(df)
}

#' Validate representations for gkcm_suffStat; assumes that df
#' has already been validated.
#'
#' @noRd
.validate_suffStat_representations <- function(representations, df) {

  if (!is.list(representations)) {
    stop(
      "'representations' must be a named list.",
      call. = FALSE
    )
  }

  column_names <- names(df)
  rep_names <- names(representations)
  n <- nrow(df)

  if (length(representations) != ncol(df) ||
      is.null(rep_names) ||
      anyNA(rep_names) ||
      any(!nzchar(rep_names)) ||
      anyDuplicated(rep_names) ||
      !setequal(rep_names, column_names)) {
    stop(
      paste0(
        "'representations' must contain exactly one named entry ",
        "for every column of 'df'."
      ),
      call. = FALSE
    )
  }

  for (column_name in column_names) {
    rep <- representations[[column_name]]

    if (!is.list(rep) ||
        length(rep) != 2 ||
        is.null(names(rep)) ||
        anyNA(names(rep)) ||
        anyDuplicated(names(rep)) ||
        !setequal(names(rep), c("type", "value"))) {
      stop(
        sprintf(
          "The representation for '%s' must contain 'type' and 'value'.",
          column_name
        ),
        call. = FALSE
      )
    }

    valid_type <- identical(rep$type, "data") ||
      identical(rep$type, "features") ||
      identical(rep$type, "gram")

    if (!valid_type) {
      stop(
        sprintf(
          paste0(
            "The representation type for '%s' must be ",
            "'data', 'features', or 'gram'."
          ),
          column_name
        ),
        call. = FALSE
      )
    }

    value <- rep$value
    missing <- is.na(df[[column_name]])
    observed <- !missing

    if (rep$type == "data") {
      if (!is.factor(value) || !is.null(dim(value))) {
        stop(
          sprintf(
            "The data representation for '%s' must be a factor vector.",
            column_name
          ),
          call. = FALSE
        )
      }

      if (length(value) != n) {
        stop(
          sprintf(
            "The factor representation for '%s' must have length %d.",
            column_name,
            n
          ),
          call. = FALSE
        )
      }

      if (any(is.na(value) != missing)) {
        stop(
          sprintf(
            paste0(
              "The missing values in the factor representation for '%s' ",
              "must match the missing values in 'df'."
            ),
            column_name
          ),
          call. = FALSE
        )
      }

      next
    }

    if (!is.matrix(value) ||
        !is.numeric(value) ||
        is.complex(value)) {
      stop(
        sprintf(
          "The representation for '%s' must be a real numeric matrix.",
          column_name
        ),
        call. = FALSE
      )
    }

    if (nrow(value) != n) {
      stop(
        sprintf(
          "The representation for '%s' must have %d rows.",
          column_name,
          n
        ),
        call. = FALSE
      )
    }

    if (rep$type == "features" && ncol(value) < 1) {
      stop(
        sprintf(
          "The feature matrix for '%s' must have at least one column.",
          column_name
        ),
        call. = FALSE
      )
    }

    if (rep$type == "gram" && ncol(value) != n) {
      stop(
        sprintf(
          "The Gram matrix for '%s' must have %d rows and columns.",
          column_name,
          n
        ),
        call. = FALSE
      )
    }

    # NaN is also considered missing by is.na(), so check it separately.
    if (any(is.nan(value)) || any(is.infinite(value))) {
      stop(
        sprintf(
          paste0(
            "The representation for '%s' must not contain ",
            "NaN or infinite values."
          ),
          column_name
        ),
        call. = FALSE
      )
    }

    if (rep$type == "features") {
      missing_rows <- value[missing, , drop = FALSE]
      observed_rows <- value[observed, , drop = FALSE]

      if (any(!is.na(missing_rows)) || anyNA(observed_rows)) {
        stop(
          sprintf(
            paste0(
              "The feature matrix for '%s' must have NA in every missing ",
              "row and finite values in every observed row."
            ),
            column_name
          ),
          call. = FALSE
        )
      }

    } else {
      missing_rows <- value[missing, , drop = FALSE]
      missing_columns <- value[, missing, drop = FALSE]
      observed_block <- value[observed, observed, drop = FALSE]

      if (any(!is.na(missing_rows)) ||
          any(!is.na(missing_columns)) ||
          anyNA(observed_block)) {
        stop(
          sprintf(
            paste0(
              "The Gram matrix for '%s' must have NA in every missing ",
              "row and column and finite values in the observed block."
            ),
            column_name
          ),
          call. = FALSE
        )
      }

      symmetric <- isTRUE(all.equal(
        observed_block,
        t(observed_block),
        tolerance = sqrt(.Machine$double.eps),
        check.attributes = FALSE
      ))

      if (!symmetric) {
        stop(
          sprintf(
            "The observed block of the Gram matrix for '%s' must be symmetric.",
            column_name
          ),
          call. = FALSE
        )
      }
    }
  }

  invisible(representations)
}

#' Validate per-column ranger arguments for gkcm_suffStat
#'
#' @noRd
.validate_suffStat_ranger_args <- function(ranger_args, df, representations) {

  if (is.null(ranger_args)) {
    ranger_args <- list()
  }

  checkmate::assert_list(
    ranger_args,
    names = if (length(ranger_args) > 0L) "named" else NULL,
    .var.name = "ranger_args"
  )
  checkmate::assert_true(!is.data.frame(ranger_args), .var.name = "ranger_args")

  if (anyDuplicated(names(ranger_args))) {
    stop("'ranger_args' cannot contain duplicate column names.", call. = FALSE)
  }

  unknown_names <- setdiff(names(ranger_args), names(df))
  if (length(unknown_names) > 0L) {
    stop(
      sprintf(
        "Unknown column names in 'ranger_args': %s.",
        paste(unknown_names, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  for (column_name in names(ranger_args)) {
    column_args <- ranger_args[[column_name]]
    if (is.null(column_args) ||
        (is.list(column_args) && !is.data.frame(column_args) && length(column_args) == 0L)) {
      stop(
        sprintf(
          paste0(
            "The `ranger_args` entry for column '%s' must be a non-empty named list; ",
            "omit this column from `ranger_args` to use the defaults."
          ),
          column_name
        ),
        call. = FALSE
      )
    }

    ranger_args[[column_name]] <- .validate_ranger_args(
      args = column_args,
      var = column_name,
      feature_type = representations[[column_name]]$type
    )
  }

  ranger_args
}

#' Validate a gkcm_suffStat object
#'
#' @noRd
.validate_gkcm_suffStat <- function(x) {

  if (!inherits(x, "gkcm_suffStat") || !is.list(x)) {
    stop(
      "'x' must be a gkcm_suffStat object.",
      call. = FALSE
    )
  }

  expected_fields <- c(
    "df",
    "representations",
    "ranger_args",
    "num_RFFs",
    "num_threads",
    "seed"
  )

  if (!identical(names(x), expected_fields)) {
    stop(
      "The gkcm_suffStat object has an invalid structure.",
      call. = FALSE
    )
  }

  .validate_suffStat_df(x$df)

  .validate_num_RFFs(x$num_RFFs)

  if (is.null(x$num_threads)) {
    stop(
      paste0(
        "The 'num_threads' field must contain a positive integer. ",
        "Create a new object with gkcm_suffStat()."
      ),
      call. = FALSE
    )
  }
  .validate_num_threads(x$num_threads)

  .validate_seed(x$seed)

  .validate_suffStat_representations(
    representations = x$representations,
    df = x$df
  )
  .validate_suffStat_ranger_args(
    ranger_args = x$ranger_args,
    df = x$df,
    representations = x$representations
  )

  invisible(x)
}

#' Only feature matrices and data are checked; not Gram matrices.
#' @noRd
.representation_is_constant <- function(rep) {

  if (rep$type == "data") {
    return(length(unique(rep$value)) == 1)
  }

  if (rep$type == "features") {
    reference <- rep$value[rep.int(1, nrow(rep$value)),,drop = FALSE]
    return(all(rep$value == reference))
  }

  FALSE
}

#' Save kernel representations for numeric data;
#' return factor for categorical data.
#' @noRd
.transform_data <- function(rep, num_RFFs) {

  if (is.numeric(rep$value)) {
    if (!is.null(num_RFFs)) {
      rep$value <- .rff(rep$value, num_RFFs)
      rep$type <- "features"
    } else {
      rep$value <- .gaussian_kernel(rep$value)
      rep$type <- "gram"
    }
  } else {
    if (!is.factor(rep$value)) {
      rep$value <- as.factor(rep$value)
    } else {
      rep$value <- droplevels(rep$value)
    }
  }

  rep
}
