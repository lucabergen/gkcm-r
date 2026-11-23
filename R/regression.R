#' Subtracts the estimated feature mean conditional on x and returns either a
#' feature or Gram matrix as the residual representation.
#' Currently we do not support linear kernels. Hence y_rep$type == "data"
#' must be a factor vector.
#'
#' @noRd
.regress_representation <- function(y_rep, x, num_threads = 1, seed, ranger_args) {

  # GKCM handles this option itself; it is not an argument of either backend.
  oob.predictions <- FALSE
  if (!is.null(ranger_args) && "oob.predictions" %in% names(ranger_args)) {
    oob.predictions <- .validate_oob_predictions(ranger_args[["oob.predictions"]])
    ranger_args[["oob.predictions"]] <- NULL
  }

  if (oob.predictions && !is.null(ranger_args) && isFALSE(ranger_args$replace) &&
      isTRUE(ranger_args$sample.fraction == 1)) {
    stop("`replace = FALSE` and `sample.fraction = 1` leave no out-of-bag observations for GKCM.",
         call. = FALSE)
  }

  if (y_rep$type == "data"){

    # ranger with gini splitrule is equivalent to full-rank kernel splitting
    # with dirac kernel, but uses category sums as sufficient statistics, which
    # can be updated faster during split search

    defaults <- list(
      mtry = ncol(x),
      num.trees = 500,
      min.node.size = 5,
      min.bucket = 1,
      max.depth = 0,
      replace = TRUE,
      sample.fraction = 1,
      respect.unordered.factors = "ignore"
    )

    args <- .merge_ranger_defaults(defaults, ranger_args)

    args <- c(args, list(
      y = y_rep$value,
      x = x,
      probability = TRUE,
      splitrule = "gini",
      seed = seed,
      num.threads = num_threads,
      oob.error = oob.predictions,
      write.forest = TRUE
    ))

    rf <- do.call(ranger::ranger, args, envir = asNamespace("ranger"))

    P <- if (oob.predictions) {
      rf$predictions
    } else {
      .predict_training_forest(rf, x, "ranger", num_threads, seed)
    }

    if (!is.matrix(P) || !all(is.finite(P))) {
      stop(if (oob.predictions) {
        "Some observations do not have finite OOB predictions. Use more trees or lower the sample fraction."
      } else {
        "The forest did not return finite training predictions."
      }, call. = FALSE)
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
      mtry = ncol(x),
      num.trees = 500,
      min.node.size = 5,
      min.bucket = 1,
      max.depth = 0,
      replace = TRUE,
      sample.fraction = 1,
      respect.unordered.factors = "ignore",
      approx.tolerance = 0.001
    )

    args <- .merge_ranger_defaults(defaults, ranger_args)

    args <- c(args, list(
      Phi = y_rep$value,
      x = x,
      num.threads = num_threads,
      seed = seed,
      oob.error = oob.predictions,
      write.forest = TRUE
    ))

    rf <- do.call(okranger::okranger, args, envir = asNamespace("okranger"))

    predictions <- if (oob.predictions) {
      rf$predictions
    } else {
      .predict_training_forest(rf, x, "okranger", num_threads, seed)
    }
    if (!is.matrix(predictions) || !identical(dim(predictions), dim(y_rep$value)) ||
        !all(is.finite(predictions))) {
      stop(if (oob.predictions) {
        "Some observations do not have finite OOB predictions. Use more trees or lower the sample fraction."
      } else {
        "The forest did not return finite training predictions with the expected dimensions."
      }, call. = FALSE)
    }
    E <- y_rep$value - predictions

    return(list(type = "features", value = E))

  } else {

    defaults <- list(
      mtry = ncol(x),
      num.trees = 500,
      min.node.size = 5,
      min.bucket = 1,
      max.depth = 0,
      replace = TRUE,
      sample.fraction = 1,
      respect.unordered.factors = "ignore",
      approx.tolerance = 0.001
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

    W <- okranger::getForestWeights(rf, training.data = x, oob = oob.predictions)
    # W may be sparse; matrix multiplication avoids densifying it.
    weight_sum <- as.vector(W %*% rep(1, ncol(W)))
    if (any(!is.finite(weight_sum) | weight_sum <= 0)) {
      stop(if (oob.predictions) {
        "Some observations do not have finite OOB weights. Use more trees or lower the sample fraction."
      } else {
        "Some observations do not have finite training-prediction weights."
      }, call. = FALSE)
    }
    B <- as.matrix(W %*% y_rep$value)
    R <- y_rep$value - B - t(B) + t(as.matrix(W %*% t(B)))

    return(list(type = "gram", value = R))
  }

}

#' Predict on the training covariates using the fitting backend's own method.
#' Both packages use the same S3 class name, so select the correct method
#' explicitly instead of relying on which package registered it last.
#' @noRd
.predict_training_forest <- function(rf, x, backend, num_threads, seed) {
  predict_forest <- utils::getFromNamespace("predict.ranger.forest", backend)
  predict_forest(rf$forest, data = x, num.threads = num_threads,
                 seed = seed)$predictions
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
