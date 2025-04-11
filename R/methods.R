#####################################################################################
#' The generalized covariance measure test of Shah and Peters.
#'
#' \code{GCM} is a function carrying out the GCM test based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param fitting_X_on_Z The fitting method for the regression X on Z.
#' @param fitting_Y_on_Z The fitting method for the regression Y on Z.
#'
#' @return A named list with fields \code{test_stat} and \code{left_side_p_value},
#' \code{right_side_p_value}, \code{both_side_p_value}.
#'
#' @examples
#' n <- 20; p <- 2
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- GCM(data, X_on_Z_fam, Y_on_Z_fam,
#'                fitting_X_on_Z = 'rf',
#'                fitting_Y_on_Z = 'glm')
#'
#' @export
GCM <- function(data, X_on_Z_fam, Y_on_Z_fam,
                fitting_X_on_Z = 'glm',
                fitting_Y_on_Z = 'glm') {

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  fitted_vals <- fit_models(data = data,
                            X_on_Z_fam = X_on_Z_fam,
                            Y_on_Z_fam = Y_on_Z_fam,
                            fitting_X_on_Z = fitting_X_on_Z,
                            fitting_Y_on_Z = fitting_Y_on_Z)

  X_on_Z_fit_vals <- fitted_vals$X_on_Z_fit_vals
  Y_on_Z_fit_vals <- fitted_vals$Y_on_Z_fit_vals

  # compute the products of residuals for each observation
  prod_resids <- (X - X_on_Z_fit_vals)*(Y - Y_on_Z_fit_vals)

  # compute the test statistic
  test_stat <- 1/sqrt(n)*sum(prod_resids)/stats::sd(prod_resids) * sqrt(n/(n-1))

  # return test statistic, GCM p-values, and related quantities
  return(list(test_stat = test_stat,
              p.left = stats::pnorm(test_stat, lower.tail = TRUE),
              p.right = stats::pnorm(test_stat, lower.tail = FALSE),
              p.both = 2*stats::pnorm(abs(test_stat), lower.tail = FALSE),
              NB.disp.param = fitted_vals$additional_info$NB.disp.param,
              unnormalized_test_stat = 1/sqrt(n)*sum(prod_resids)))
}


#####################################################################################
#' The distilled conditional randomization test.
#'
#' \code{dCRT} is a function carrying out the dCRT based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param fitting_X_on_Z The fitting method for the regression X on Z.
#' @param fitting_Y_on_Z The fitting method for the regression Y on Z.
#' @param B The number of resamples to draw.
#'
#' @return A named list with fields \code{test_stat}, \code{left_side_p_value},
#' \code{right_side_p_value}, and \code{both_side_p_value}.
#'
#' @examples
#' n <- 80; p <- 2; B <- 100
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- dCRT(data, X_on_Z_fam, Y_on_Z_fam,
#'                 fitting_X_on_Z = 'rf',
#'                 fitting_Y_on_Z = 'glm',
#'                 B = 2000)
#'
#' @export
dCRT <- function(data, X_on_Z_fam, Y_on_Z_fam,
                 fitting_X_on_Z = 'glm',
                 fitting_Y_on_Z = 'glm',
                 B = 2000) {

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  fitted_vals <- fit_models(data = data,
                            X_on_Z_fam = X_on_Z_fam,
                            Y_on_Z_fam = Y_on_Z_fam,
                            fitting_X_on_Z = fitting_X_on_Z,
                            fitting_Y_on_Z = fitting_Y_on_Z)

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

  # return test statistic, dCRT p-values, and related quantities
  return(list(test_stat = test_stat,
              p.left = p.left,
              p.right = p.right,
              p.both = 2*min(p.left, p.right),
              NB.disp.param = fitted_vals$additional_info$NB.disp.param))
}


#####################################################################################
#' The saddlepoint approximation to the dCRT.
#'
#' \code{spaCRT} is a function carrying out the saddlepoint approximation to the
#' dCRT based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param fitting_X_on_Z The fitting method for the regression X on Z.
#' @param fitting_Y_on_Z The fitting method for the regression Y on Z.
#' @param R Upper bound of search space for the saddlepoint
#'
#' @return A named list with fields \code{test_stat}, \code{left_side_p_value},
#' \code{right_side_p_value}, \code{both_side_p_value}, and
#' \code{gcm.default}.
#' gcm.default returns TRUE if GCM was employed due to the failure of spaCRT.
#'
#' @examples
#' n <- 50; p <- 4
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rbinom(n = n, size = 1, prob = 0.7),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "binomial"
#' spaCRT(data, X_on_Z_fam, Y_on_Z_fam,
#'        fitting_X_on_Z = 'rf',
#'        fitting_Y_on_Z = 'glm')
#'
#' @export
spaCRT <- function(data, X_on_Z_fam, Y_on_Z_fam,
                   fitting_X_on_Z = 'glm',
                   fitting_Y_on_Z = 'glm',
                   R = 5) {

  fitted_vals <- fit_models(data = data,
                            X_on_Z_fam = X_on_Z_fam,
                            Y_on_Z_fam = Y_on_Z_fam,
                            fitting_X_on_Z = fitting_X_on_Z,
                            fitting_Y_on_Z = fitting_Y_on_Z) |> suppressWarnings()

  spa_result <- spa_cdf(X = data$X, Y = data$Y,
                        X_on_Z_fit_vals = fitted_vals$X_on_Z_fit_vals,
                        Y_on_Z_fit_vals = fitted_vals$Y_on_Z_fit_vals,
                        fam = X_on_Z_fam,
                        R = abs(R),
                        max_expansions = 10) |> suppressWarnings()

  NB.disp.param <- fitted_vals$additional_info$NB.disp.param

  return(spa_result |> append(list(NB.disp.param = NB.disp.param), after = 4))
}







# spacrt - methods.R

