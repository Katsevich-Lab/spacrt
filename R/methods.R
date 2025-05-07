#####################################################################################
#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM} is a function carrying out the GCM test based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param X
#'   Numeric vector of length n for the predictor variable.
#' @param Y
#'   Numeric vector of length n for the response variable.
#' @param Z
#'   Numeric matrix (n × p) of covariates.
#' @param family
#'   Named list with elements `XZ` and `YZ` specifying the distribution or loss
#'   (e.g. `"binomial"`, `"poisson"`) for each model.
#'   Ignored for any model where you supply your own fitted values via `fitted.own`.
#' @param method
#'   Named list with elements `XZ` and `YZ` that selects the fitting engine
#'   for each model. Each element must be a string (e.g. `"glm"`, `"rf"`).
#'   Ignored for any model where you supply your own fitted values via `fitted.own`.
#' @param fitted.own
#'   Named list with elements `XZ` and `YZ` of user-supplied fitted values
#'   (numeric vectors of length n). For each non-NULL element, that model is
#'   treated as custom and neither `method` nor `family` is used.
#' @param alternative
#'   A character string specifying the alternative hypothesis for the test.
#'   Values can be "two.sided" (default), "greater", or "less".
#'
#' @return
#' A named list containing the following fields:
#' \describe{
#'   \item{test_stat}{The test statistic.}
#'   \item{p_value}{The p-value under the specified alternative.}
#' }
#'
#' @examples
#' n <- 20; p <- 2
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' GCM(data, X_on_Z_fam, Y_on_Z_fam,
#'     fitting_X_on_Z = 'rf',
#'     fitting_Y_on_Z = 'glm')
#'
#' @export
GCM <- function(X, Y, Z,
                family,
                method = list(XZ = 'glm', YZ = 'glm'),
                fitted.own = list(XZ = NULL, YZ = NULL),
                alternative = 'two.sided') {

  n <- length(X)

  fitted_vals <- fit_models(X, Y, Z,
                            family = list(XZ = family$XZ,
                                          YZ = family$YZ),
                            method = list(XZ = method$XZ,
                                          YZ = method$YZ),
                            fitted.own = list(XZ = fitted.own$XZ,
                                              YZ = fitted.own$XZ)) |> suppressWarnings()

  X_on_Z_fit_vals <- fitted_vals$X_on_Z_fit_vals
  Y_on_Z_fit_vals <- fitted_vals$Y_on_Z_fit_vals

  # compute the products of residuals for each observation
  prod_resids <- (X - X_on_Z_fit_vals)*(Y - Y_on_Z_fit_vals)

  # compute the test statistic
  test_stat <- 1/sqrt(n)*sum(prod_resids)/stats::sd(prod_resids) * sqrt(n/(n-1))

  # compute p-value based on sideness
  p.left <- stats::pnorm(test_stat, lower.tail = TRUE)
  p.right <- stats::pnorm(test_stat, lower.tail = FALSE)
  p.both <- 2 * stats::pnorm(abs(test_stat), lower.tail = FALSE)

  pval <- switch(alternative,
                 less = c(p.left = p.left),
                 greater = c(p.right = p.right),
                 two.sided = c(p.both = p.both),
                 stop("Invalid value for `alternative`"))

  # return test statistic and GCM p-value
  return(list(test_stat = test_stat,
              p_value = pval))
}


#####################################################################################
#' The distilled conditional randomization test.
#'
#' \code{dCRT} is a function carrying out the dCRT based on GLMs for `X|Z` and `Y|Z`.
#'
#' @inheritParams GCM
#' @param B The number of resamples to draw (Default value is 2000).
#'
#' @return
#' An named list containing the following fields:
#' \describe{
#'   \item{test_stat}{The test statistic.}
#'   \item{p_value}{The p-value under the specified alternative.}
#' }
#'
#' @examples
#' n <- 80; p <- 2; B <- 100
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#'
#' dCRT(data, X_on_Z_fam, Y_on_Z_fam,
#'      fitting_X_on_Z = 'rf',
#'      fitting_Y_on_Z = 'glm',
#'      B = 2000)
#'
#' @export
dCRT <- function(X, Y, Z,
                 family,
                 method = list(XZ = 'glm', YZ = 'glm'),
                 fitted.own = list(XZ = NULL, YZ = NULL),
                 alternative = 'two.sided',
                 B = 2000) {

  n <- length(X)

  fitted_vals <- fit_models(X, Y, Z,
                            family = list(XZ = family$XZ,
                                          YZ = family$YZ),
                            method = list(XZ = method$XZ,
                                          YZ = method$YZ),
                            fitted.own = list(XZ = fitted.own$XZ,
                                              YZ = fitted.own$XZ)) |> suppressWarnings()

  X_on_Z_fit_vals <- fitted_vals$X_on_Z_fit_vals
  Y_on_Z_fit_vals <- fitted_vals$Y_on_Z_fit_vals

  # compute the products of residuals for each observation
  prod_resids <- (X - X_on_Z_fit_vals)*(Y - Y_on_Z_fit_vals)
  # compute the test statistic
  test_stat <- 1/sqrt(n) * sum(prod_resids)

  prod_resid_resamp <- c()

  for(b in 1:B){
    # resampling X from X|Z
    resamp_X <- dCRT_dist(n = n,
                          fitted.val = X_on_Z_fit_vals,
                          fam = family$XZ)

    # compute the products of residuals for each resampled observation
    prod_resid_resamp[b] <- 1/sqrt(n) * sum((resamp_X - X_on_Z_fit_vals)*(Y - Y_on_Z_fit_vals))
  }

  # compute p-values
  p.left <- 1/(B+1) * (1 + sum(prod_resid_resamp <= test_stat))
  p.right <- 1/(B+1) * (1 + sum(prod_resid_resamp >= test_stat))
  p.both <- 2 * min(c(p.left, p.right))

  pval <- switch(alternative,
                 less = c(p.left = p.left),
                 greater = c(p.right = p.right),
                 two.sided = c(p.both = p.both),
                 stop("Invalid value for `alternative`"))

  # return test statistic and dCRT p-value
  return(list(test_stat = test_stat,
              p_value = pval))
}


#####################################################################################
#' The saddlepoint approximation to the dCRT.
#'
#' \code{spaCRT} is a function carrying out the saddlepoint approximation to the
#' dCRT based on GLMs for `X|Z` and `Y|Z`.
#'
#' @inheritParams GCM
#'
#' @return
#' A named list containing the following fields:
#' \describe{
#'   \item{\code{test_stat}}{The test statistic.}
#'   \item{\code{p_value}}{The p-value under the specified alternative.}
#'   \item{\code{spa.success}}{A logical variable that returns TRUE if the
#'   saddlepoint equation could be solved; otherwise, the backup method (GCM)
#'   was employed due to the failure of spaCRT.}
#' }
#'
#' @examples
#' ## Example 1
#' n <- 50; p <- 4
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#'
#' spaCRT(data, X_on_Z_fam, Y_on_Z_fam,
#'        fitting_X_on_Z = 'rf',
#'        fitting_Y_on_Z = 'glm',
#'        alternative = 'greater')
#'
#' ## Example 2
#' n <- 100; p <- 10
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.7),
#'              Y = rbinom(n = n, size = 1, prob = 0.2),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "binomial"
#'
#' dtrain <- xgboost::xgb.DMatrix(data = data$Z, label = data$Y)
#' model.Y <- xgboost::xgboost(data = dtrain,
#'                             objective = "binary:logistic",
#'                             nrounds = 50, verbose = 0)
#' predicted <- stats::predict(model.Y, newdata = data$Z)
#'
#' spaCRT(data, X_on_Z_fam, Y_on_Z_fam,
#'        fitting_X_on_Z = 'glm',
#'        fitting_Y_on_Z = 'own',
#'        fit_vals_Y_on_Z_own = predicted)
#'
#' @export
spaCRT <- function(X, Y, Z,
                   family,
                   method = list(XZ = 'glm', YZ = 'glm'),
                   fitted.own = list(XZ = NULL, YZ = NULL),
                   alternative = 'two.sided') {

  fitted_vals <- fit_models(X, Y, Z,
                            family = list(XZ = family$XZ,
                                          YZ = family$YZ),
                            method = list(XZ = method$XZ,
                                          YZ = method$YZ),
                            fitted.own = list(XZ = fitted.own$XZ,
                                              YZ = fitted.own$XZ)) |> suppressWarnings()

  spa_result <- spa_cdf(X = X, Y = Y,
                        X_on_Z_fit_vals = fitted_vals$X_on_Z_fit_vals,
                        Y_on_Z_fit_vals = fitted_vals$Y_on_Z_fit_vals,
                        fam = family$XZ,
                        R = 5) |> suppressWarnings()

  # compute p-value based on alternative
  p.left <- spa_result$p.left
  p.right <- spa_result$p.right
  p.both <- spa_result$p.both

  pval <- switch(alternative,
                 less = c(p.left = p.left),
                 greater = c(p.right = p.right),
                 two.sided = c(p.both = p.both),
                 stop("Invalid value for `alternative`"))

  # return test statistic and spaCRT p-value
  return(list(test_stat = spa_result$test_stat,
              p_value = pval,
              spa.success = spa_result$spa.success))
}







# spacrt - methods.R

