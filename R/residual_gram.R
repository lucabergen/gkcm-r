
#' @noRd
.regress <- function(y, X, feat_repr_y = NULL) {

  if (is.factor(y)){

    # Regress the explicit feature map of categorical targets

    # Use in-bag predictions
    # rf <- ranger::ranger(formula  = y ~ ., data = cbind(y = y, X = X),
    #                      mtry = ncol(X), num.trees = ncol(X)*100,
    #                      probability = TRUE, splitrule = "gini")
    # pred <- predict(rf, data = data.frame(X = X), type = "response")
    # P <- pred$predictions

    # Use out-of-bag predictions
    rf <- ranger::ranger(formula  = y ~ ., data = cbind(y = y, X = X),
                         mtry = floor(ncol(X)/2), num.trees = ncol(X)*50,
                         probability = TRUE, splitrule = "gini", keep.inbag = F)
    P <- rf$predictions # n x k matrix

    # Ensure label levels match columns of P
    classes <- colnames(P)
    y <- factor(y, levels = classes)

    # Probability residuals
    R <- feat_repr_y - P

    # Inner products
    R %*% t(R)

  } else {

    # Implicitly form and regress the feature map of numeric targets
    rf <- drf::drf(X, y, num.trees = ncol(X)*100, bandwidth = 1,
                   honesty = F, response.scaling = F,
                   min.node.size = 5, mtry = 100)
    # Use newdata = Z as additional argument to get in-sample (non-OOB) weights
    W <- drf::get_sample_weights(rf) # |> as.matrix()

    .compute_drf_residuals(feat_repr_y, W)

  }

}

#' @noRd
.get_residual_gram <- function(y, X, feat_repr_y = NULL) {

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
      .regress(y, X, .get_feat_repr(y))
    } else {
      .regress(y, X, feat_repr_y)
    }

  }
}

#' @noRd
.regress <- function(y, X, feat_repr_y = NULL) {

  if (is.factor(y)){

    # Regress the explicit feature map of categorical targets

    # Use in-bag predictions
    # rf <- ranger::ranger(formula  = y ~ ., data = cbind(y = y, X = X),
    #                      mtry = ncol(X), num.trees = ncol(X)*100,
    #                      probability = TRUE, splitrule = "gini")
    # pred <- predict(rf, data = data.frame(X = X), type = "response")
    # P <- pred$predictions

    # Use out-of-bag predictions
    rf <- ranger::ranger(formula  = y ~ ., data = cbind(y = y, X = X),
                         mtry = ncol(X), num.trees = ncol(X)*100,
                         probability = TRUE, splitrule = "gini", keep.inbag = F)
    P <- rf$predictions # n x k matrix

    # Ensure label levels match columns of P
    classes <- colnames(P)
    y <- factor(y, levels = classes)

    # Probability residuals
    R <- feat_repr_y - P

    # Inner products
    R %*% t(R)

  } else {

    # Implicitly form and regress the feature map of numeric targets
    rf <- drf::drf(X, y, num.trees = ncol(X)*100, bandwidth = 1,
                   honesty = F, response.scaling = F,
                   min.node.size = 5, mtry = 100)
    # Use newdata = Z as additional argument to get in-sample (non-OOB) weights
    W <- drf::get_sample_weights(rf) # |> as.matrix()

    .compute_drf_residuals(feat_repr_y, W)

  }

}

.get_residual_gram <- function(y, X, feat_repr_y = NULL) {

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
      .regress(y, X, .get_feat_repr(y))
    } else {
      .regress(y, X, feat_repr_y)
    }

  }
}


