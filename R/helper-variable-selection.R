#' Compute all conditional means for Y|X_{-j}
#'
#' @param fitted_model Fitted regression model
#' @param x Query point
#' @param conditional_prob Conditional probability matrix (p-by-M)
#' @param support_x Support set of X
#' @param lambda Lambda regularization choice in cv.glmnet
#' @param post_lasso A logical value. If TRUE, post_lasso should be used and the
#' `fitted_model` should be a low-dimensional GLM
#'
#' @return Conditional probability vector E(Y|X_{-j} = x_{-j})
#' @export
compute_all_means <- function(fitted_model, x, conditional_prob, support_x, lambda = "min",
                              post_lasso = FALSE){

  # extract the dimension of the problem
  p <- length(x)

  # obtain a square matrix
  x_impute <- matrix(rep(x, p), nrow = p, ncol = p, byrow = TRUE)

  # separate the analysis for different post_lasso parameters
  if(post_lasso){

    # compute the all the conditional expectations Y|X_{-j}
    fix_one_fitted <- sapply(sort(support_x), function(s){

      # obtain the fitted model on (x_1,..,x_{j-1}, \{0,\ldots,M\},x_{j+1},...,x_p)
      diag(x_impute) <- s

      # separate the case when act_set is NULL
      if(length(fitted_model$act_set) == 0){
        # use just intercept term
        stats::predict(fitted_model$model,
                       newdata = data.frame(rep(1, nrow(x_impute))),
                       type = "response")
      }else{
        # use predict function to obtain the leave-one-out fitted values
        stats::predict(fitted_model$model,
                       newdata = data.frame(x_impute)[, sort(fitted_model$act_set), drop = FALSE],
                       type = "response")
      }

    })

  }else{

    # compute the all the conditional expectations Y|X_{-j}
    fix_one_fitted <- sapply(sort(support_x), function(s){

      # obtain the fitted model on (x_1,..,x_{j-1}, \{0,\ldots,M\},x_{j+1},...,x_p)
      diag(x_impute) <- s

      # use predict function to obtain the leave-one-out fitted values
      stats::predict(fitted_model, newx = x_impute, s = sprintf("lambda.%s", lambda),
                     type = "response")

    })
  }

  # compute the
  integrate_one_fitted <- rowSums(as.matrix(fix_one_fitted) * conditional_prob)

  # return the output
  return(integrate_one_fitted)
}


#' Compute the conditional mean using power trick (efficient version)
#'
#' @inheritParams compute_all_means
#' @param X data matrix of dim n-by-p
#' @param conditional_prob_mat Conditional probability matrix of dim n-by-(p*M) and M is the size of support_x
#'
#' @return Matrix of dim n-by-p
#' @export
compute_all_means_efficient <- function(fitted_model, X, conditional_prob_mat,
                                        support_x, lambda = "min"){

  # extract the dimension of the problem
  p <- ncol(X)
  n <- nrow(X)

  # predicted matrix to be imputed
  predicted_mat <- matrix(0, nrow = n, ncol = p)

  # extract the nonzero component in the model
  act_set <- sort(which(as.vector(stats::coef(fitted_model, s = sprintf("lambda.%s", lambda)))[-1] != 0))
  predicted_mat[, dplyr::setdiff(1:p, act_set)] <- stats::predict(fitted_model,
                                                                  s = sprintf("lambda.%s", lambda),
                                                                  newx = X, type = "response")

  # separate the case when act_set is of length zero or not
  if (length(act_set) > 0){

    # loop over the act_set
    integrate_one_fitted <- sapply(act_set, function(act_beta){

      # compute the all the conditional expectations Y|X_{-j} for j in act_set
      fix_one_fitted <- sapply(sort(support_x), function(s){

        # impute the X to be the value in support_x
        x_impute <- X
        x_impute[, act_beta] <- s

        # use predict function to obtain the leave-one-out fitted values
        stats::predict(fitted_model, newx = x_impute, s = sprintf("lambda.%s", lambda),
                       type = "response")

      })

      # compute the integrated prediction
      rowSums(fix_one_fitted * conditional_prob_mat[, act_beta + p * (0 : (length(support_x) - 1))])
    })

    # finish the imputation
    predicted_mat[, act_set] <- integrate_one_fitted
  }

  # return the output
  return(predicted_mat)
}


#' This is a post-lasso fitting function
#'
#' @param X A matrix including n rows and p columns
#' @param Y A vector of length n
#' @param family A GLM family
#' @param lambda_param_vec A vector including either "lambda.min" or "lambda.1se"
#'
#' @return A list of fitted models and active set selected from lasso algorithm
#' @export
post_lasso <- function(X, Y, family = "binomial",
                       lambda_param_vec = c("lambda.min", "lambda.1se")){

  # fit Y on X using lasso
  lasso_model <- glmnet::cv.glmnet(x = X, y = Y, family = family)

  # transform the data to data frame
  dimnames(X) <- list(
    sample = 1:nrow(X),
    predictor = 1:ncol(X)
  )

  # extract the lambda.min and lambda.1se parameters
  fitted_model <- list(
    lambda.1se = list(NULL),
    lambda.min = list(NULL)
  )
  for (lambda_param in lambda_param_vec) {

    # extract the non-zero coefficient
    act_set <- which(as.vector(stats::coef(lasso_model, s = lambda_param))[-1] != 0)
    fitted_model[[lambda_param]]$act_set <- act_set

    # perform the glm.fit
    if(length(act_set) == 0){
      glm_fitted <- stats::glm(Y ~ 1, data = data.frame(Y),
                               family = family)
    }else{
      X_act <- X[, sort(act_set), drop = FALSE]
      glm_fitted <- stats::glm(Y ~ ., data = data.frame(Y, X_act),
                               family = family)
    }

    # extract the fitted model
    fitted_model[[lambda_param]]$model <- glm_fitted
  }

  # return the output
  return(fitted_model)
}
