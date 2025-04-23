#####################################################################################
#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM} is a function carrying out the GCM test based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data A (non-empty) named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z} (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param fitting_X_on_Z The fitting method for the regression X on Z.
#' (values can be \code{glm} (default), \code{rf}, \code{prob_forest}, or \code{own})
#' @param fitting_Y_on_Z The fitting method for the regression Y on Z.
#' (values can be \code{glm} (default), \code{rf}, \code{prob_forest}, or \code{own})
#' @param fit_vals_X_on_Z_own Vector of fitted values for X on Z in case the user's custom method.
#' Works only if fitting_X_on_Z = 'own'.
#' @param fit_vals_Y_on_Z_own Vector of fitted values for Y on Z in case the user's custom method.
#' Works only if fitting_Y_on_Z = 'own'.
#' @param alternative A character string specifying the alternative hypothesis
#' (values can be "two.sided" (default), "greater", or "less"; meant for both-sided, right-sided,
#' and left-sided p-values, respectively).
#'
#' @return A named list with fields \code{test_stat} and \code{p_value}.
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
GCM <- function(data, X_on_Z_fam, Y_on_Z_fam,
                fitting_X_on_Z = 'glm',
                fitting_Y_on_Z = 'glm',
                fit_vals_X_on_Z_own = NULL,
                fit_vals_Y_on_Z_own = NULL,
                alternative = 'two.sided') {

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  fitted_vals <- fit_models(data = data,
                            X_on_Z_fam = X_on_Z_fam,
                            Y_on_Z_fam = Y_on_Z_fam,
                            fitting_X_on_Z = fitting_X_on_Z,
                            fitting_Y_on_Z = fitting_Y_on_Z,
                            fit_vals_X_on_Z_own = fit_vals_X_on_Z_own,
                            fit_vals_Y_on_Z_own = fit_vals_Y_on_Z_own) |> suppressWarnings()

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
#' @return A named list with fields \code{test_stat} and \code{p_value}.
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
dCRT <- function(data, X_on_Z_fam, Y_on_Z_fam,
                 fitting_X_on_Z = 'glm',
                 fitting_Y_on_Z = 'glm',
                 fit_vals_X_on_Z_own = NULL,
                 fit_vals_Y_on_Z_own = NULL,
                 alternative = 'two.sided',
                 B = 2000) {

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  fitted_vals <- fit_models(data = data,
                            X_on_Z_fam = X_on_Z_fam,
                            Y_on_Z_fam = Y_on_Z_fam,
                            fitting_X_on_Z = fitting_X_on_Z,
                            fitting_Y_on_Z = fitting_Y_on_Z,
                            fit_vals_X_on_Z_own = fit_vals_X_on_Z_own,
                            fit_vals_Y_on_Z_own = fit_vals_Y_on_Z_own) |> suppressWarnings()

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
                          fam = X_on_Z_fam)

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
#' @return A named list with fields \code{test_stat}, \code{p_value}, and \code{spa.success}.
#' \code{spa.success} returns TRUE if the saddlepoint equation could be solved; otherwise,
#' the backup method (GCM) was employed.
#'
#' @examples
#' n <- 50; p <- 4
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rbinom(n = n, size = 1, prob = 0.7),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "poisson"
#' Y_on_Z_fam <- "binomial"
#' spaCRT(data, X_on_Z_fam, Y_on_Z_fam,
#'        fitting_X_on_Z = 'glm',
#'        fitting_Y_on_Z = 'glm',
#'        alternative = 'greater')
#'
#' n <- 100; p <- 10
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.3),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#'
#' pois_obj <- function(b) sum(exp(data$Z %*% b) - data$Y * (data$Z %*% b))
#' M <- matrix(0, nrow = ncol(data$Z), ncol = 1)
#' beta_hat <- stats::optim(M, pois_obj, method = "BFGS")$par
#' fit_vals_Y_on_Z_own <- as.numeric(exp(data$Z %*% beta_hat))
#'
#' spaCRT(data, X_on_Z_fam, Y_on_Z_fam,
#'        fitting_X_on_Z = 'rf',
#'        fitting_Y_on_Z = 'own',
#'        fit_vals_Y_on_Z_own = fit_vals_Y_on_Z_own)
#'
#' @export
spaCRT <- function(data, X_on_Z_fam, Y_on_Z_fam,
                   fitting_X_on_Z = 'glm',
                   fitting_Y_on_Z = 'glm',
                   fit_vals_X_on_Z_own = NULL,
                   fit_vals_Y_on_Z_own = NULL,
                   alternative = 'two.sided') {

  fitted_vals <- fit_models(data = data,
                            X_on_Z_fam = X_on_Z_fam,
                            Y_on_Z_fam = Y_on_Z_fam,
                            fitting_X_on_Z = fitting_X_on_Z,
                            fitting_Y_on_Z = fitting_Y_on_Z,
                            fit_vals_X_on_Z_own = fit_vals_X_on_Z_own,
                            fit_vals_Y_on_Z_own = fit_vals_Y_on_Z_own) |> suppressWarnings()

  spa_result <- spa_cdf(X = data$X, Y = data$Y,
                        X_on_Z_fit_vals = fitted_vals$X_on_Z_fit_vals,
                        Y_on_Z_fit_vals = fitted_vals$Y_on_Z_fit_vals,
                        fam = X_on_Z_fam,
                        R = 5, max_expansions = 10) |> suppressWarnings()

  NB.disp.param <- fitted_vals$additional_info$NB.disp.param

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

