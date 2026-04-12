
#' @noRd
.check_var <- function(x, n = NULL) {

  if (is.null(n)) {
    if (length(x) < 10) {
      return(sprintf("must have at least 10 observations (has %d)", length(x)))
    }
  } else {
    if (length(x) != n) {
      return("number of observations doesn't match Z")
    }
  }

  # must be a 1D atomic vector
  if (!isTRUE(checkmate::check_atomic_vector(x, any.missing = FALSE))) {
    return("must be a 1D atomic vector with no missing values")
  }

  # numeric: additionally forbid Inf/-Inf/NaN
  if (checkmate::test_numeric(x)) {
    res <- checkmate::check_numeric(x, any.missing = FALSE, finite = TRUE)
    if (!isTRUE(res)) return(res)
    if (length(unique(x)) == 1) return("must not be constant")
    return(TRUE)
  }

  # other allowed atomic types (any.missing already enforced above)
  if (checkmate::test_integer(x))   { if (length(unique(x)) == 1) return("must not be constant"); return(TRUE) }
  if (checkmate::test_logical(x))   { if (length(unique(x)) == 1) return("must not be constant"); return(TRUE) }
  if (checkmate::test_factor(x))    { if (length(unique(as.integer(x))) == 1) return("must not be constant"); return(TRUE) }
  if (checkmate::test_character(x)) { if (length(unique(x)) == 1) return("must not be constant"); return(TRUE) }

  "must be a vector (complex and raw not allowed)"

}

#' @noRd
.assert_var <- checkmate::makeAssertionFunction(.check_var)

.check_cond <- function(Z) {

  # NULL is allowed (marginal test)
  if (is.null(Z)) return(TRUE)

  # Enforce minimum of 10 observations
  if (is.data.frame(Z) || is.matrix(Z)) {
    n <- nrow(Z)
  } else {
    n <- length(Z)
  }

  if (n < 10) {
    return(sprintf("must have at least 10 observations (has %d)", n))
  }

  # data.frame (mixed types)
  if (is.data.frame(Z)) {
    res <- checkmate::check_data_frame(
      Z,
      types = c("numeric", "integer", "logical", "factor", "character"),
      any.missing = FALSE
    )
    if (!isTRUE(res)) return(res)

    # numeric columns must be finite (no Inf/-Inf/NaN)
    num_cols <- vapply(Z, is.numeric, logical(1))
    if (any(num_cols) && any(!is.finite(as.matrix(Z[num_cols])))) {
      return("numeric columns must contain only finite values (no Inf/-Inf/NaN)")
    }

    # columns must not be constant
    const_col <- vapply(Z, function(col) length(unique(col)) == 1, logical(1))
    if (any(const_col)) {
      return(sprintf("columns must not be constant (constant: %s)", names(Z)[which(const_col)[1]]))
    }

    return(TRUE)
  }

  # matrix (numeric or integer)
  if (is.matrix(Z)) {

    if (is.integer(Z)) {
      res <- checkmate::check_matrix(Z, mode = "integer", any.missing = FALSE)
      if (!isTRUE(res)) return(res)
      if (any(apply(Z, 2, function(col) length(unique(col)) == 1))) {
        return("columns must not be constant")
      }
      return(TRUE)
    }

    if (is.numeric(Z)) {
      if (any(!is.finite(Z))) {
        return("numeric matrix must contain only finite values (no Inf/-Inf/NaN)")
      }
      res <- checkmate::check_matrix(Z, mode = "numeric", any.missing = FALSE)
      if (!isTRUE(res)) return(res)
      if (any(apply(Z, 2, function(col) length(unique(col)) == 1))) {
        return("columns must not be constant")
      }
      return(TRUE)
    }

    return("matrix must be numeric or integer")
  }

  # vector (atomic and allowed types only)
  if (!isTRUE(checkmate::check_atomic_vector(Z, any.missing = FALSE))) {
    return("must be NULL or a 1D vector, numeric/integer matrix, or data.frame")
  }

  if (checkmate::test_numeric(Z)) {
    res <- checkmate::check_numeric(Z, any.missing = FALSE, finite = TRUE)
    if (!isTRUE(res)) return(res)
    if (length(unique(Z)) == 1) return("must not be constant")
    return(TRUE)
  }

  if (checkmate::test_integer(Z))   { if (length(unique(Z)) == 1) return("must not be constant"); return(TRUE) }
  if (checkmate::test_logical(Z))   { if (length(unique(Z)) == 1) return("must not be constant"); return(TRUE) }
  if (checkmate::test_factor(Z))    { if (length(unique(as.integer(Z))) == 1) return("must not be constant"); return(TRUE) }
  if (checkmate::test_character(Z)) { if (length(unique(Z)) == 1) return("must not be constant"); return(TRUE) }

  "must be NULL (marginal test) or a vector, matrix, or data.frame (no complex or raw values)"

}

#' @noRd
.assert_cond <- checkmate::makeAssertionFunction(.check_cond)

#' @noRd
.check_n <- function(x, y) {
  nx <- length(x)
  ny <- length(y)

  if (nx != ny) {
    return(sprintf("x and y must have the same number of observations (x has %d, y has %d)", nx, ny))
  }

  TRUE
}

#' @noRd
.assert_n <- checkmate::makeAssertionFunction(.check_n)

#' @noRd
.check_df <- function(df) {

  # Must be a data.frame
  if (!is.data.frame(df)) {
    return("must be a data.frame")
  }

  # Minimum number of observations
  n <- nrow(df)
  if (n < 10) {
    return(sprintf("must have at least 10 observations (has %d)", n))
  }

  # Minimum number of columns
  p <- ncol(df)
  if (p < 2) {
    return(sprintf("must have at least 2 columns (has %d)", p))
  }

  # Only allow these column types, and no missing values
  res <- checkmate::check_data_frame(
    df,
    types = c("numeric", "integer", "logical", "factor", "character"),
    any.missing = FALSE
  )
  if (!isTRUE(res)) return(res)

  # Numeric columns must be finite (no Inf/-Inf/NaN)
  num_cols <- vapply(df, is.numeric, logical(1))
  if (any(num_cols) && any(!is.finite(as.matrix(df[num_cols])))) {
    return("numeric columns must contain only finite values (no Inf/-Inf/NaN)")
  }

  # columns must not be constant
  const_col <- vapply(df, function(col) length(unique(col)) == 1, logical(1))
  if (any(const_col)) {
    return(sprintf("columns must not be constant (constant: %s)", names(df)[which(const_col)[1]]))
  }

  TRUE

}

#' @noRd
.assert_df <- checkmate::makeAssertionFunction(.check_df)

#' @noRd
.check_num_threads <- function(num_threads) {

  if (is.null(num_threads)) {
    return(TRUE)
  }

  # Must be a single positive integer
  if (!is.numeric(num_threads) || length(num_threads) != 1 || is.na(num_threads)) {
    return("must be NULL or a single positive integer")
  }

  if (!is.finite(num_threads) || num_threads < 1 || num_threads != as.integer(num_threads)) {
    return("must be NULL or a single positive integer")
  }

  # Must not exceed available cores
  n_avail <- parallelly::availableCores()
  if (num_threads > n_avail) {
    return(sprintf("must not exceed available cores (%d)", n_avail))
  }

  TRUE
}

#' @noRd
.assert_num_threads <- checkmate::makeAssertionFunction(.check_num_threads)
