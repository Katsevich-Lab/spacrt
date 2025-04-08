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
#' @param fitting_method The fitting method for the regressions X on Z and Y on Z.
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
#'                fitting_method = 'glm')
#'
#' @export
GCM <- function(data, X_on_Z_fam, Y_on_Z_fam,
                fitting_method = 'glm') {

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  fitted_vals <- fit_models(data = data,
                            fitting_method = fitting_method,
                            X_on_Z_fam = X_on_Z_fam, Y_on_Z_fam = Y_on_Z_fam)

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
#' @param fitting_method The fitting method for the regressions X on Z and Y on Z.
#' @param B The number of resamples to draw.
#' @param normalize A logical variable indicating whether the dCRT should be
#' normalized.
#' @param return_resamples A logical variable indicating whether to return the
#' resampled test statistics.
#'
#' @return A named list with fields \code{test_stat}, \code{left_side_p_value},
#' \code{right_side_p_value}, \code{both_side_p_value} and
#' \code{resamples}. Here, \code{resamples} is a vector of length \code{B}. It
#' is returned only if \code{return_resamples == TRUE}.
#'
#' @examples
#' n <- 80; p <- 2; B <- 100; normalize <- FALSE; return_resamples <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' results <- dCRT(data, X_on_Z_fam, Y_on_Z_fam,
#'                 fitting_method = 'glm', B = 2000,
#'                 normalize, return_resamples)
#'
#' @export
dCRT <- function(data, X_on_Z_fam, Y_on_Z_fam,
                 fitting_method = 'glm',
                 B = 2000, normalize = FALSE, return_resamples = FALSE) {

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  fitted_vals <- fit_models(data = data,
                            fitting_method = fitting_method,
                            X_on_Z_fam = X_on_Z_fam, Y_on_Z_fam = Y_on_Z_fam)

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
#' @param fitting_method The fitting method for the regressions X on Z and Y on Z.
#' @param normalize A logical variable indicating whether the spaCRT is based on
#' the normalized test statistic.
#' @param return_cdf A logical variable indicating whether to return the CDF
#' @param R Upper bound of search space for the saddlepoint
#'
#' @return A named list with fields \code{test_stat}, \code{left_side_p_value},
#' \code{right_side_p_value}, \code{both_side_p_value}, \code{cdf}, and
#' \code{gcm.default}.
#' Here, cdf is a function that takes in a value t and returns the
#' saddlepoint approximation to the CDF of the resampling distribution of the
#' test statistic evaluated at t. This function is returned only if
#' return_cdf == TRUE.
#' gcm.default returns TRUE if GCM was employed due to the failure of spaCRT.
#'
#' @examples
#' n <- 50; p <- 6; normalize <- FALSE; return_cdf <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "poisson"
#' spaCRT(data, X_on_Z_fam, Y_on_Z_fam, fitting_method = 'glm', normalize)
#'
#' n <- 50; p <- 4; normalize <- FALSE; return_cdf <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rbinom(n = n, size = 1, prob = 0.7),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "binomial"
#' spaCRT(data, X_on_Z_fam, Y_on_Z_fam, fitting_method = "rf")
#'
#' @export
spaCRT <- function(data, X_on_Z_fam, Y_on_Z_fam,
                   fitting_method = 'glm',
                   normalize = FALSE, return_cdf = FALSE,
                   R = 5) {

  fitted_vals <- fit_models(data = data,
                            X_on_Z_fam = X_on_Z_fam, Y_on_Z_fam = Y_on_Z_fam,
                            fitting_method = fitting_method) |> suppressWarnings()

  spa_result <- spa_cdf(X = data$X, Y = data$Y,
                        X_on_Z_fit_vals = fitted_vals$X_on_Z_fit_vals,
                        Y_on_Z_fit_vals = fitted_vals$Y_on_Z_fit_vals,
                        fam = X_on_Z_fam,
                        R = abs(R),
                        max_expansions = 10) |> suppressWarnings()

  NB.disp.param <- fitted_vals$additional_info$NB.disp.param

  return(spa_result |> append(list(NB.disp.param = NB.disp.param), after = 4))
}




#####################################################################################
#' Score Test
#'
#' \code{score.test} is a function carrying out the saddlepoint approximation to the
#' dCRT based on GLMs for `X|Z` and `Y|Z`.
#'
#' @param data A named list with fields \code{X} (an nx1 vector for the predictor
#' variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
#' (an nxp matrix of covariates).
#' @param X_on_Z_fam The GLM family for the regression of X on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param Y_on_Z_fam The GLM family for the regression of Y on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#'
#' @examples
#' n <- 50; p <- 2; normalize <- FALSE; return_cdf <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X_on_Z_fam <- "binomial"
#' Y_on_Z_fam <- "negative.binomial"
#' results <- score.test(data, X_on_Z_fam, Y_on_Z_fam)
#' results$test_stat
#'
#' @return A named list with fields \code{test_stat} and \code{left_side_p_value},
#' \code{right_side_p_value} and \code{both_side_p_value}.
#'
#' @export
score.test <- function(data, X_on_Z_fam, Y_on_Z_fam){

  # extract (X,Y,Z) from inputted data
  X <- data$X; Y <- data$Y; Z <- data$Z
  n <- length(X)

  # fit X on Z regression
  X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = X_on_Z_fam))

  # fit Y on Z regression
  if(Y_on_Z_fam == "negative.binomial"){
    # First try to fit the model using glm.nb
    temp.result <- tryCatch({
      Y_on_Z_fit <- suppressWarnings(MASS::glm.nb(Y ~ Z))
      NB.disp.param <- Y_on_Z_fit$theta
      list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
    },
    error = function(e) {
      aux_info_Y_on_Z <- nb_precomp(list(Y = Y, Z = Z))

      Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z,
                              family = MASS::negative.binomial(aux_info_Y_on_Z$theta_hat),
                              mustart = aux_info_Y_on_Z$fitted_values))
      NB.disp.param <- aux_info_Y_on_Z$theta_hat

      list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
    })
  }else if(Y_on_Z_fam == 'poisson'){
    aux_info_Y_on_Z <- nb_precomp(list(Y = Y, Z = Z))

    Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z,
                                              family = stats::poisson(),
                                              mustart = aux_info_Y_on_Z$fitted_values))
    NB.disp.param <- NA

    temp.result <- list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
  }else{
    Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = Y_on_Z_fam))
    NB.disp.param <- NA

    temp.result <- list(Y_on_Z_fit = Y_on_Z_fit, NB.disp.param = NB.disp.param)
  }

  # perform score test
  test_stat <- statmod::glm.scoretest(fit = temp.result$Y_on_Z_fit, x2 = X)

  # return test statistic, score test p-values, and related quantities
  return(list(test_stat = test_stat,
              p.left = stats::pnorm(test_stat, lower.tail = TRUE),
              p.right = stats::pnorm(test_stat, lower.tail = FALSE),
              p.both = 2*stats::pnorm(abs(test_stat), lower.tail = FALSE),
              NB.disp.param = temp.result$NB.disp.param))
}




# spacrt - methods.R

