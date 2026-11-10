
#' Computes test statistic and eigenvalues of T matrix.
#' @noRd
.compute_test_components <- function(eps_rep, xi_rep) {

  n <- nrow(eps_rep$value)

  if (eps_rep$type != "gram" && xi_rep$type != "gram") {
    # Compute components using feature matrices
    V <- .rowwise_outer_product(eps_rep$value, xi_rep$value)
    C <- crossprod(eps_rep$value, xi_rep$value)
    statistic <- sum(C^2) / n

    M <- (tcrossprod(V) - tcrossprod(rowSums(V)) / n) / (n - 1)
    M <- (M + t(M)) / 2
    eigenvalues <- eigen(M, symmetric = TRUE, only.values = TRUE)$values

  } else {
    # Ensure both residual representations are Gram matrices
    if (eps_rep$type != "gram") {
      eps_rep$value <- tcrossprod(eps_rep$value)
    } else if (xi_rep$type != "gram") {
      xi_rep$value <- tcrossprod(xi_rep$value)
    }

    # Compute components using Gram representations
    R <- eps_rep$value * xi_rep$value
    statistic <- sum(R) / n
    M <- .center_gram(R) / (n - 1)
    eigenvalues <- eigen(M, symmetric = TRUE, only.values = TRUE)$values

  }

  return(list(statistic = statistic, eigenvalues = eigenvalues))

}

#' @noRd
.compute_p_value <- function(components) {

  # Filter nonzero eigenvalues
  tol <- 100 * .Machine$double.eps * max(1, max(abs(components$eigenvalues)))

  eigenvalues<- components$eigenvalues
  eigenvalues[eigenvalues < 0 & abs(eigenvalues) <= tol] <- 0
  eigenvalues <- eigenvalues[eigenvalues > tol]

  if (length(eigenvalues) == 0) {
    p_val <- ifelse(components$statistic <= tol, 1, 0)
  } else if (components$statistic <= 0) {
    p_val <- 1
  } else {
    dvs <- CompQuadForm::davies(components$statistic, lambda = eigenvalues)
    if (dvs$ifault == 0 && is.finite(dvs$Qq) && dvs$Qq >= 0 && dvs$Qq <= 1) {
      p_val <- dvs$Qq
    } else {
      imh <- CompQuadForm::imhof(components$statistic, lambda = eigenvalues)
      p_val <- min(1, max(0, imh$Qq))
    }
  }

  return(p_val)
}

#' @noRd
.rowwise_outer_product <- function(R,S) {

  p <- ncol(R)
  q <- ncol(S)
  n <- nrow(R)

  V <- vapply(seq_len(n),
              function(i) {
                as.vector(R[i,, drop = FALSE] %o% S[i,, drop = FALSE])
              },
              numeric(p*q))

  matrix(V, nrow = p*q, ncol = n)

}
