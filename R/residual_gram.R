
#' @noRd
.regress <- function(y, X, feat_repr_y = NULL, num_threads = 1) {

  if (is.factor(y)){

    # Regress the explicit feature map of categorical targets

    # Use in-bag predictions
    rf <- ranger::ranger(formula  = y ~ ., data = cbind(y = y, X = X),
                         mtry = ncol(X), num.trees = ncol(X)*50,
                         probability = TRUE, splitrule = "gini",
                         num.threads = num_threads)
    pred <- stats::predict(rf, data = data.frame(X = X),
                    type = "response",
                    num.threads = num_threads)
    P <- pred$predictions # n x k matrix

    # Ensure label levels match columns of P
    classes <- colnames(P)
    # Reorder and restrict to match ranger output (unused levels are dropped)
    feat_repr_y <- feat_repr_y[, classes, drop = FALSE]

    # Probability residuals
    R <- feat_repr_y - P

    # Inner products
    R %*% t(R)

  } else {

    # Implicitly form and regress the feature map of numeric targets
    rf <- drf::drf(X = X, Y = y, num.trees = ncol(X)*100, bandwidth = 1,
                   honesty = F, response.scaling = F,
                   min.node.size = 5, mtry = 100,
                   num.threads = num_threads, compute.oob.predictions = F)

    # Use newdata = X to get in-sample (non-OOB) weights
    W <- stats::predict(rf, newdata = X, num.threads = num_threads)$weights

    .compute_drf_residuals(feat_repr_y, W)

  }

}

#' @noRd
.get_residual_gram <- function(y, X, feat_repr_y = NULL, num_threads = 1) {

  if (is.null(X)) {

    # Center (subtract marginal mean of implicit features)
    if (is.null(feat_repr_y) || is.factor(y)) {
      L <- .get_gram(y)
    } else {
      L <- feat_repr_y
    }
    .center_gram(L)

  } else {

    # Regress (subtract conditional mean of implicit or explicit features)
    if (is.null(feat_repr_y)) {
      .regress(y, X, .get_feat_repr(y), num_threads)
    } else {
      .regress(y, X, feat_repr_y, num_threads)
    }

  }
}
