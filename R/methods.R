#####################################################################################
#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM} is a function carrying out the GCM test based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param X
#'   Numeric vector of length \eqn{n}, representing the predictor variable.
#' @param Y
#'   Numeric vector of length \eqn{n}, representing the response variable.
#' @param Z
#'   Numeric matrix with \eqn{n} rows and \eqn{p} columns, representing covariates.
#' @param family
#'   Named list with elements `XZ` and `YZ` specifying the model family for \eqn{X \mid Z} and
#'   \eqn{Y \mid Z} for each model. Each list element must be a string (e.g. `"binomial"`,
#'   `"poisson"`). Ignored for any model where you supply your own fitted values via `fitted`.
#' @param method
#'   Named list with elements `XZ` and `YZ` that selects the model-fitting method to use
#'   for each model. Each element must be a string (e.g. `"glm"`, `"random_forest"`).
#'   Ignored for any model where you supply your own fitted values via `fitted`.
#' @param fitted
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
#' @seealso \code{\link{supported_models}}
#' @examples
#' n <- 200; p <- 4
#' set.seed(1234)
#' X <- rbinom(n = n, size = 1, prob = 0.3)
#' Y <- rpois(n = n, lambda = 1)
#' Z <- matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p)
#'
#' # fit both models via GLM
#' res.GCM.1 <- GCM(X = X, Y = Y, Z = Z,
#'                  family = list(XZ = "binomial", YZ = "poisson"),
#'                  method = list(XZ = "glm", YZ = "glm"))
#' res.GCM.1
#'
#' # custom fit for Y|Z, explicit method/family only for X|Z
#' user_fit_Y <- glm(Y ~ Z,
#'                   family = poisson(),
#'                   data = data.frame(Y = Y, Z = Z))$fitted.values |> unname()
#'
#' res.GCM.2 <- GCM(X = X, Y = Y, Z = Z,
#'                  method = list(XZ = "random_forest"),
#'                  family = list(XZ = "binomial"),
#'                  fitted = list(XZ = NULL, YZ = user_fit_Y),
#'                  alternative = "greater")
#' res.GCM.2
#'
#' @export
GCM <- function(X, Y, Z,
                family,
                method,
                fitted = list(XZ = NULL, YZ = NULL),
                alternative = 'two.sided') {

  check_inputs_main(X, Y, Z,
                    family, method, fitted,
                    alternative,
                    func = 'GCM')

  n <- length(X)

  fitted_vals <- fit_models(X, Y, Z,
                            family = list(XZ = family$XZ,
                                          YZ = family$YZ),
                            method = list(XZ = method$XZ,
                                          YZ = method$YZ),
                            fitted = list(XZ = fitted$XZ,
                                              YZ = fitted$YZ)) |> suppressWarnings()

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
#' @param B The number of resamples to draw (Default value is 5000).
#'
#' @return
#' An named list containing the following fields:
#' \describe{
#'   \item{test_stat}{The test statistic.}
#'   \item{p_value}{The p-value under the specified alternative.}
#' }
#' @seealso \code{\link{supported_models}}
#' @examples
#' n <- 200; p <- 4
#' set.seed(1234)
#' X <- rbinom(n = n, size = 1, prob = 0.3)
#' Y <- rpois(n = n, lambda = 1)
#' Z <- matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p)
#'
#' # fit both models via GLM
#' res.dCRT.1 <- dCRT(X = X, Y = Y, Z = Z,
#'                    family = list(XZ = "binomial", YZ = "poisson"),
#'                    method = list(XZ = "glm", YZ = "glm"))
#' res.dCRT.1
#'
#' # custom fit for Y|Z, explicit method/family only for X|Z
#' user_fit_Y <- glm(Y ~ Z,
#'                   family = poisson(),
#'                   data = data.frame(Y = Y, Z = Z))$fitted.values |> unname()
#'
#' res.dCRT.2 <- dCRT(X = X, Y = Y, Z = Z,
#'                    method = list(XZ = "random_forest"),
#'                    family = list(XZ = "binomial"),
#'                    fitted = list(XZ = NULL, YZ = user_fit_Y),
#'                    alternative = "greater")
#' res.dCRT.2
#'
#' @export
dCRT <- function(X, Y, Z,
                 family,
                 method = list(XZ = 'glm', YZ = 'glm'),
                 fitted = list(XZ = NULL, YZ = NULL),
                 alternative = 'two.sided',
                 B = 5000) {

  check_inputs_main(X, Y, Z,
                    family, method, fitted,
                    alternative,
                    func = 'dCRT', B)

  n <- length(X)

  fitted_vals <- fit_models(X, Y, Z,
                            family = list(XZ = family$XZ,
                                          YZ = family$YZ),
                            method = list(XZ = method$XZ,
                                          YZ = method$YZ),
                            fitted = list(XZ = fitted$XZ,
                                              YZ = fitted$YZ)) |> suppressWarnings()

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
#'   \item{test_stat}{The test statistic.}
#'   \item{p_value}{The p-value under the specified alternative.}
#'   \item{spa.success}{A logical variable that returns TRUE if the
#'   saddlepoint equation could be solved; otherwise, the backup method (GCM)
#'   was employed due to the failure of spaCRT.}
#' }
#' @seealso \code{\link{supported_models}}
#' @examples
#' n <- 200; p <- 4
#' set.seed(1234)
#' X <- rbinom(n = n, size = 1, prob = 0.3)
#' Y <- rpois(n = n, lambda = 1)
#' Z <- matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p)
#'
#' # fit both models via GLM
#' res.spaCRT.1 <- spaCRT(X = X, Y = Y, Z = Z,
#'                        family = list(XZ = "binomial", YZ = "poisson"),
#'                        method = list(XZ = "glm", YZ = "glm"))
#' res.spaCRT.1
#'
#' # custom fit for Y|Z, explicit method/family only for X|Z
#' dtrain <- xgboost::xgb.DMatrix(data = Z, label = Y)
#' model.Y <- xgboost::xgboost(data = dtrain,
#'                             objective = "count:poisson",
#'                             nrounds = 100, verbose = 0)
#' user_fit_Y <- stats::predict(model.Y, newdata = Z)
#'
#' res.spaCRT.2 <- spaCRT(X = X, Y = Y, Z = Z,
#'                        method = list(XZ = "random_forest"),
#'                        family = list(XZ = "binomial"),
#'                        fitted = list(XZ = NULL, YZ = user_fit_Y),
#'                        alternative = "greater")
#' res.spaCRT.2
#'
#' @export
spaCRT <- function(X, Y, Z,
                   family,
                   method = list(XZ = 'glm', YZ = 'glm'),
                   fitted = list(XZ = NULL, YZ = NULL),
                   alternative = 'two.sided') {

  check_inputs_main(X, Y, Z,
                    family, method, fitted,
                    alternative,
                    func = 'spaCRT')

  fitted_vals <- fit_models(X, Y, Z,
                            family = list(XZ = family$XZ,
                                          YZ = family$YZ),
                            method = list(XZ = method$XZ,
                                          YZ = method$YZ),
                            fitted = list(XZ = fitted$XZ,
                                              YZ = fitted$YZ)) |> suppressWarnings()

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
