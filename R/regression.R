#' Subtracts the estimated feature mean conditional on x and returns either a
#' feature or Gram matrix as the residual representation.
#' Currently we do not support linear kernels. Hence y_rep$type == "data"
#' must be a factor vector.
#'
#' @noRd
.regress_representation <- function(y_rep, x, num_threads = 1, seed, ranger_args) {

  if (y_rep$type == "data"){

    # ranger with gini splitrule is equivalent to full-rank kernel splitting
    # with dirac kernel, but uses category sums as sufficient statistics, which
    # can be updated faster during split search

    defaults <- list(
      mtry = ncol(x),
      num.trees = 200,
      min.node.size = 10,
      max.depth = 0
    )

    args <- .merge_ranger_defaults(defaults, ranger_args)

    args <- c(args, list(
      y = y_rep$value,
      x = x,
      probability = TRUE,
      splitrule = "gini",
      seed = seed,
      num.threads = num_threads,
      oob.error = TRUE
    ))

    rf <- do.call(ranger::ranger, args, envir = asNamespace("ranger"))

    P <- rf$predictions

    if (!all(is.finite(P))) {
      stop("Some observations do not have finite OOB predictions. Use more trees.")
    }

    class_idx <- match(as.character(y_rep$value), colnames(P))

    if (anyNA(class_idx)) {
      stop("At least one class is not predicted by ranger.")
    }

    n <- length(y_rep$value)
    Y <- matrix(0, nrow = n, ncol = ncol(P))
    Y[cbind(seq_len(n), class_idx)] <- 1

    E <- Y - P

    return(list(type = "features", value = E))

  } else if (y_rep$type == "features") {

    defaults <- list(
      approx.tolerance = 0.001,
      mtry = ncol(x),
      num.trees = 200,
      min.node.size = 5,
      max.depth = 0
    )

    args <- .merge_ranger_defaults(defaults, ranger_args)

    args <- c(args, list(
      Phi = y_rep$value,
      x = x,
      num.threads = num_threads,
      seed = seed,
      oob.error = TRUE
    ))

    rf <- do.call(okranger::okranger, args, envir = asNamespace("okranger"))

    E <- y_rep$value - rf$predictions

    return(list(type = "features", value = E))

  } else {

    defaults <- list(
      approx.tolerance = 0.001,
      mtry = ncol(x),
      num.trees = 200,
      min.node.size = 5,
      max.depth = 0
    )

    args <- .merge_ranger_defaults(defaults, ranger_args)

    args <- c(args, list(
      K = y_rep$value,
      x = x,
      num.threads = num_threads,
      seed = seed,
      keep.inbag = TRUE,
      oob.error = FALSE,
      write.forest = TRUE
    ))

    rf <- do.call(okranger::okranger, args, envir = asNamespace("okranger"))

    W <- okranger::getForestWeights(rf, training.data = x, oob = TRUE)
    B <- as.matrix(W %*% y_rep$value)
    R <- y_rep$value - B - t(B) + t(as.matrix(W %*% t(B)))

    return(list(type = "gram", value = R))
  }

}

#' Use defaults for tuning parameters not overridden by ranger_args
#'
#' @noRd
.merge_ranger_defaults <- function(defaults, ranger_args) {

  if (is.null(ranger_args)) {
    return(defaults)
  }

  defaults[names(ranger_args)] <- ranger_args
  defaults
}

#' Used in gkcm_indepTest. Derives a deterministic regression seed from the
#' global seed given to gkcm_indepTest, the index integers, and the target
#' string (must be either "x" or "y").
#'
#' @noRd
.derive_regression_seed <- function(seed, x, y, S, target) {

  key <- paste(
    seed,
    paste(sort(c(x, y)), collapse = ","),
    paste(sort(S), collapse = ","),
    target,
    sep = "|"
  )

  hash <- digest::digest2int(key)

  # Handle the integer value represented as NA in R
  if (is.na(hash)) {
    hash <- 0
  }

  # Map the result to a positive integer
  as.integer(as.double(hash) %% .Machine$integer.max + 1)
}
