
#' Chooses and labels representation and signals choice to user when ambiguous.
#' Features are preferred to Gram and Gram is preferred to data.
#' @noRd
.choose_representation <- function(args, var) {
  if (var == "x") {
    priority <- c("x", "K", "Phi")
  } else if (var == "y") {
    priority <- c("y", "L", "Psi")
  } else {
    stop("`var` must be either 'x' or 'y'.")
  }

  if (!is.list(args) || is.null(names(args))) {
    stop("`args` must be a named list.")
  }

  if (!all(priority %in% names(args))) {
    stop(
      sprintf(
        "`args` must include: %s.",
        paste(sprintf("'%s'", priority), collapse = ", ")
      )
    )
  }

  supplied <- priority[!vapply(args[priority], is.null, logical(1))]

  if (length(supplied) == 0L) {
    stop(
      sprintf(
        "No representation was supplied for variable %s. Supply one of: %s.",
        var,
        paste(sprintf("'%s'", priority), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  selected <- supplied[length(supplied)]

  if (length(supplied) > 1L) {
    warning(
      sprintf(
        "Multiple representations were supplied for variable %s via arguments %s; using '%s'.",
        var,
        paste(sprintf("'%s'", supplied), collapse = ", "),
        selected
      ),
      call. = FALSE
    )
  }

  if (selected %in% c("x", "y")) {
    type <- "data"
  } else if (selected %in% c("Phi", "Psi")) {
    type <- "features"
  } else {
    type <- "gram"
  }

  list(type = type, value = args[[selected]])
}

#' @noRd
.rff <- function(x, n_features = 20, lengthscale = 1) {

  x <- as.numeric(x)
  m <- as.integer(n_features / 2L)

  omega <- stats::rnorm(m, mean = 0, sd = 1 / lengthscale)
  projection <- outer(x, omega)

  sqrt(1 / m) * cbind(cos(projection), sin(projection))
}

#' @noRd
.onehot_encode <- function(x){

  classes <- levels(x)
  n <- length(x)
  k <- length(classes)

  OH <- matrix(0L, nrow = n, ncol = k)
  colnames(OH) <- classes

  OH[cbind(seq_len(n), as.integer(x))] <- 1L

  OH

}

#' @noRd
.gaussian_kernel <- function(x, sigma = 1) {

  # For vector-valued x consider parallelDist::parDist
  x <- as.vector(x)
  exp(-outer(x, x, "-")^2 / (2 * sigma^2))

}

#' Subtracts marginal feature mean and returns either a feature or Gram matrix
#' as the residual representation.
#' Currently we do not support linear kernels. Hence y_rep$type == "data"
#' must be a factor vector.
#'
#' @noRd
.center_representation <- function(y_rep) {

  n <- NROW(y_rep$value)

  if (y_rep$type == "data") {

    class_idx <- as.integer(y_rep$value)
    k <- nlevels(y_rep$value)
    p <- tabulate(class_idx, nbins = k) / n
    E <- matrix(-p, nrow = n, ncol = k, byrow = TRUE)
    E[seq_len(n) + (class_idx - 1) * n] <- 1 - p[class_idx]
    colnames(E) <- levels(y_rep$value)

    return(list(type = "features", value = E))

  } else if (y_rep$type == "features") {

    E <- sweep(y_rep$value, 2, colMeans(y_rep$value), "-")

    return(list(type = "features", value = E))

  } else {

    R <- .center_gram(y_rep$value)

    return(list(type = "gram", value = R))

  }

}

#' @noRd
.center_gram <- function(K) {

  r <- rowMeans(K)
  m <- mean(K)

  # Center: Kc = K - r*1' - 1*r' + m
  Kc <- sweep(K, 1, r, FUN = "-")
  Kc <- sweep(Kc, 2, r, FUN = "-")
  Kc <- Kc + m

  # Symmetrize
  (Kc + t(Kc)) / 2

}
