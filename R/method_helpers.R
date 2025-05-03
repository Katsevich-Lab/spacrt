#' @useDynLib spacrt, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


############################################################################################
#' \code{spa_cdf} SPA to CDF of T_n = S_n / sqrt(n)
#'
#' @param X The point where the CGF will be computed.
#' @param Y The point where the CGF will be computed.
#' @param X_on_Z_fit_vals X_on_Z_fit$fitted.values
#' @param Y_on_Z_fit_vals Y_on_Z_fit$fitted.values
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param R stats::uniroot() search space endpoint
#' should be broadened.
#'
#' @return \code{test statistic}, \code{left-sided p-value}, \code{right-sided p-value},
#' \code{both-sided p-value}, and \code{spa.success} which specifies whether the saddlepoint
#' equation could be solved. If not, a backup method (GCM) had to be employed as a backup.
#'
#' @examples
#' n <- 100; p <- 2; normalize <- FALSE; return_cdf <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X <- data$X; Y <- data$Y; Z <- data$Z
#' X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = "binomial"))
#' Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = "poisson"))
#' spacrt:::spa_cdf(X = X, Y = Y,
#'                  X_on_Z_fit_vals = X_on_Z_fit$fitted.values,
#'                  Y_on_Z_fit_vals = Y_on_Z_fit$fitted.values,
#'                  fam = "binomial", R = 1000)
spa_cdf <- function(X, Y,
                    X_on_Z_fit_vals,
                    Y_on_Z_fit_vals,
                    fam,
                    R = 5){

  P <- X_on_Z_fit_vals
  W <- Y - Y_on_Z_fit_vals
  n <- length(P)

  # compute the products of residuals for each observation
  prod_resids <- (X - P) * W

  # compute the test statistic
  test_stat <- 1/sqrt(n) * sum(prod_resids)

  t <- test_stat + 1/sqrt(n) * sum(P*W)

  # current_lower <- -abs(R)
  # current_upper <- abs(R)
  # success_uniroot <- FALSE
  #
  # for (i in seq_len(max_expansions)) {
  #   tryCatch({
  #     # solve the saddlepoint equation
  #     s.hat <- stats::uniroot(
  #       f = function(s) d1_wcgf(s, P = P, W = W, fam) - sqrt(n)*t,
  #       lower = current_lower, upper = current_upper, tol = .Machine$double.eps)$root
  #
  #     success_uniroot <- TRUE
  #     break
  #   }, error = function(e) {
  #     # expand the search space if the saddlepoint is not found
  #     expansion_factor <- ifelse(i <= max_expansions/2, 2, 10)
  #     current_lower <<- current_lower * expansion_factor
  #     current_upper <<- current_upper * expansion_factor
  #   }
  #   )
  #
  #   if(success_uniroot == TRUE) break
  # }

  success_uniroot <- FALSE

  tryCatch({
    # try to solve the saddlepoint equation
    s.hat <- stats::uniroot(
      f = function(s) d1_wcgf(s, P = P, W = W, fam) - sqrt(n)*t,
      lower = current_lower, upper = current_upper,
      extendInt = "yes",
      tol = .Machine$double.eps)$root

    success_uniroot <- TRUE
  }, error = function(e) {message("stats::uniroot() failed: ", conditionMessage(e))}
  )


  if(success_uniroot == TRUE && {
      suppressWarnings({
        r.hat <- sign(s.hat) * sqrt(2 * (sqrt(n)* s.hat * t -
                                         wcgf(s = s.hat, P = P, W = W, fam)))

        # Lugannani-Rice formula
        p.left <- stats::pnorm(r.hat) + stats::dnorm(r.hat) *
          (1/r.hat - 1/(s.hat*sqrt(d2_wcgf(s = s.hat, P = P, W = W, fam))))
      })

      # decide if p.left is NA or beyond the range [0, 1] or not
      all(p.left >= 0, p.left <= 1, !is.na(p.left))
    }
    ){
       res <- list(test_stat = t - 1/sqrt(n) * sum(P*W),
                   p.left = p.left,
                   p.right = 1 - p.left,
                   p.both = 2*min(c(p.left, 1 - p.left)),
                   spa.success = TRUE)
  }else {
    test_stat <- sum(prod_resids)/(stats::sd(prod_resids) * sqrt(n-1))

    res <- list(test_stat = test_stat,
                p.left = stats::pnorm(test_stat, lower.tail = TRUE),
                p.right = stats::pnorm(test_stat, lower.tail = FALSE),
                p.both = 2*stats::pnorm(abs(test_stat), lower.tail = FALSE),
                spa.success = FALSE)
  }

  return(res)
}



############################################################################################
#' \code{dCRT_dist} is a function that returns simulated data from an appropriate
#' distribution depending on a specified  GLM family
#'
#' @param n The point where the CGF will be computed.
#' @param fitted.val A vector containing the fitted parameter values by
#' fitting a GLM to X on Z.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#'
#' @return Simulated data from an appropriate distribution.
dCRT_dist <- function(n, fitted.val, fam){

  if(fam == 'binomial') return(stats::rbinom(n = n, size = 1, prob = fitted.val))
  if(fam == 'gaussian') return()
  if(fam == 'Gamma') return()
  if(fam == 'inverse.gaussian') return()
  if(fam == 'poisson') return(stats::rpois(n = n, lambda = fitted.val))
  if(fam == 'quasi') return()
  if(fam == 'quasibinomial') return()
  if(fam == 'quasipoisson') return()
}


############################################################################################
#' \code{wcgf} is a function computing the cumulant generating function (CGF) of
#' distributions, multiplied by a weight function, from the GLM family
#'
#' @param s The point where the CGF will be computed.
#' @param P A vector containing the parameter values of the family of distributions.
#' @param W A vector containing the weights.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#'
#' @return CGF of the weighted distribution evaluated at \code{s}.
wcgf <- function(s, P, W, fam){

  if(fam == 'binomial') return(sum(log(exp(s*W)*P + 1 - P)))
  if(fam == 'gaussian') return()
  if(fam == 'Gamma') return()
  if(fam == 'inverse.gaussian') return()
  if(fam == 'poisson') return(sum(P*(exp(s*W) - 1)))
  if(fam == 'quasi') return()
  if(fam == 'quasibinomial') return()
  if(fam == 'quasipoisson') return()
}


############################################################################################
#' \code{d1_wcgf} is a function computing the derivative of the weighted cumulant
#' generating function (WCGF) of distributions from GLM family
#'
#' @inheritParams wcgf
#'
#' @return First derivative of CGF of the weighted distribution evaluated at \code{s}.
d1_wcgf <- function(s, P, W, fam){

  # if(fam == 'binomial') return(sum((W*P*exp(s*W)) / (exp(s*W)*P + 1 - P)))
  if(fam == 'binomial') return(sum((W*P) / (P + (1 - P) * exp(-s*W))))
  if(fam == 'gaussian') return()
  if(fam == 'Gamma') return()
  if(fam == 'inverse.gaussian') return()
  if(fam == 'poisson') return(sum(P * exp(s*W) * W))
  if(fam == 'quasi') return()
  if(fam == 'quasibinomial') return()
  if(fam == 'quasipoisson') return()
}


############################################################################################
#' \code{d2_wcgf} is a function computing the hessian of the weighted cumulant
#' generating function (WCGF) of distributions from GLM family
#'
#' @inheritParams wcgf
#'
#' @return Second derivative of CGF of the weighted distribution evaluated at \code{s}.
d2_wcgf <- function(s, P, W, fam){

  if(fam == 'binomial'){
    Q <- 1 - P
    # return(sum((W^2*P*Q*exp(s*W)) / (exp(s*W)*P + Q)^2))
    return(sum((W^2*P*Q) / (exp(s*W)*P^2 + 2 * P * Q + Q^2 * exp(-s*W))))
  }

  if(fam == 'gaussian'){
    Q <- 1 - P
    return()
  }

  if(fam == 'Gamma'){
    Q <- 1 - P
    return()
  }

  if(fam == 'inverse.gaussian'){
    Q <- 1 - P
    return()
  }

  if(fam == 'poisson'){
    return(sum(P*exp(s*W)*W^2))
  }

  if(fam == 'quasi'){
    Q <- 1 - P
    return()
  }

  if(fam == 'quasibinomial'){
    Q <- 1 - P
    return()
  }

  if(fam == 'quasipoisson'){
    Q <- 1 - P
    return()
  }
}


############################################################################################
#' \code{nb_precomp} is a function computing the dispersion parameter in negative
#' binomial regression
#'
#' @param V A list containing the response V
#' @param Z A list containing the covariate Z
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{fitted_values}{The fitted values from the Poisson regression of \code{Y} on \code{Z}.}
#'   \item{theta_hat}{The estimated dispersion parameter (theta) for the negative binomial model, computed via maximum likelihood or method of moments.}
#' }
nb_precomp <- function(V,Z){

  # Y <- data$Y; Z <- data$Z

  # Fit Poisson GLM: Y ~ Z
  pois_fit <- stats::glm.fit(y = V, x = Z, family = stats::poisson())
  # df <- data.frame(V = V, Z = Z)
  # pois_fit <- stats::glm(V ~ Z, data = df, family = stats::poisson())

  # Estimate NB dispersion parameter using C++ function estimate_theta
  theta_hat <- estimate_theta(
    y = Y,
    mu = pois_fit$fitted.values,
    dfr = pois_fit$df.residual,
    limit = 50,
    eps = (.Machine$double.eps)^(1/4)
  )[[1]]

  return(list(fitted_values = pois_fit$fitted.values,
              theta_hat = theta_hat))
}







############################################################################################
#' \code{fit_models} is a function carrying out the regressions `X` on `Z` and `Y` on `Z`.
#'
#' @inheritParams GCM
#'
#' @return A named list of fitted values of X|Z and Y|Z.
fit_models <- function(data,
                       X_on_Z_fam, Y_on_Z_fam,
                       fitting_X_on_Z = 'glm',
                       fitting_Y_on_Z = 'glm',
                       fit_vals_X_on_Z_own = NULL,
                       fit_vals_Y_on_Z_own = NULL){

   # extract (X,Y,Z) from inputted data
   X <- data$X; Y <- data$Y; Z <- data$Z

   # fit X on Z regression
   X_on_Z_fit_vals <- fit_single_model(V = X, Z = Z,
                                       V_on_Z_fam = X_on_Z_fam,
                                       fitting_V_on_Z = fitting_X_on_Z,
                                       fit_vals_V_on_Z_own = fit_vals_X_on_Z_own)

   # fit Y on Z regression
   Y_on_Z_fit_vals <- fit_single_model(V = Y, Z = Z,
                                       V_on_Z_fam = Y_on_Z_fam,
                                       fitting_V_on_Z = fitting_Y_on_Z,
                                       fit_vals_V_on_Z_own = fit_vals_Y_on_Z_own)

   return(list(X_on_Z_fit_vals = X_on_Z_fit_vals,
               Y_on_Z_fit_vals = Y_on_Z_fit_vals))
}



############################################################################################
#' \code{fit_single_model} is a function carrying out a single regression `V` on `Z`, where
#' `V` is the variable/response and `Z` is the covariate.
#'
#' @param V Vector of variable/response.
#' @param Z Matrix of covariates.
#' @param V_on_Z_fam The GLM family for the regression of V on Z
#' (values can be \code{gaussian}, \code{binomial}, \code{poisson}, \code{negative.binomial}, etc).
#' @param fitting_V_on_Z The fitting method for the regression V on Z
#' (values can be \code{glm} (default), \code{rf}, \code{prob_forest}, or \code{own}).
#' @param fit_vals_V_on_Z_own Vector of fitted values for V on Z in case the user's custom method.
#' Works only if fitting_V_on_Z = 'own'.
#'
#' @return A vector of fitted values of V|Z.
fit_single_model <- function(V, Z,
                             V_on_Z_fam,
                             fitting_V_on_Z = 'glm',
                             fit_vals_V_on_Z_own = NULL){

  if(fitting_V_on_Z == 'glm'){
    # fit V on Z regression when fitting method is glm
    if(V_on_Z_fam == "negative.binomial"){
      aux_info_V_on_Z <- nb_precomp(V = V, Z = Z)

      V_on_Z_fit <- stats::glm(V ~ Z,
                               family = MASS::negative.binomial(aux_info_V_on_Z$theta_hat),
                               mustart = aux_info_V_on_Z$fitted_values) |> suppressWarnings()

      V_on_Z_fit_vals <- V_on_Z_fit$fitted.values
    } else{
      V_on_Z_fit <- suppressWarnings(stats::glm(V ~ Z, family = V_on_Z_fam))
      V_on_Z_fit_vals <- V_on_Z_fit$fitted.values
    }
  } else if(fitting_V_on_Z %in% c('rf','prob_forest')){
    # fit V on Z regression when fitting method is random forest
    if(V_on_Z_fam == "binomial"){
      p.forest.V <- grf::probability_forest(X = as.matrix(Z), Y = as.factor(V))
      p.hat.V <- stats::predict(p.forest.V, as.matrix(Z), estimate.variance = F)

      V_on_Z_fit_vals <- p.hat.V$predictions[ ,"1"]
    }
  } else if(fitting_V_on_Z == 'own') {
    # Validate that fit_vals_V_on_Z_own is provided and contains necessary components
    if(!is.numeric(fit_vals_V_on_Z_own)) {
      stop("fit_vals_V_on_Z_own must be a vector containing 'V_on_Z_fit_vals' when using fitting_V_on_Z = 'own'")
    }

    # Extract pre-computed values from fit_vals_V_on_Z_own
    V_on_Z_fit_vals <- fit_vals_V_on_Z_own
  }

  return(V_on_Z_fit_vals)
}




# spacrt - method_helpers.R
