#####################################################################################
#' \code{spa_cdf} SPA to CDF of T_n = S_n / sqrt(n)
#'
#' @param X The point where the CGF will be computed.
#' @param Y The point where the CGF will be computed.
#' @param X_on_Z_fit_vals X_on_Z_fit$fitted.values
#' @param Y_on_Z_fit_vals Y_on_Z_fit$fitted.values
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param R stats::uniroot() search space endpoint
#' @param max_expansions Maximum number of times stats::uniroot() search space should be broadened
#' @return Simulated data from an appropriate distribution.
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
#' @noRd
spa_cdf <- function(X, Y, X_on_Z_fit_vals, Y_on_Z_fit_vals, fam, R, max_expansions = 10){

  P <- X_on_Z_fit_vals
  W <- Y - Y_on_Z_fit_vals
  n <- length(P)

  # compute the products of residuals for each observation
  prod_resids <- (X - P) * W

  # compute the test statistic
  test_stat <- 1/sqrt(n) * sum(prod_resids)

  t <- test_stat + 1/sqrt(n) * sum(P*W)

  current_lower <- -abs(R)
  current_upper <- abs(R)
  success_uniroot <- FALSE

  for (i in seq_len(max_expansions)) {
    tryCatch({
      # solve the saddlepoint equation
      s.hat <- stats::uniroot(
        f = function(s) d1_wcgf(s, P = P, W = W, fam) - sqrt(n)*t,
        lower = current_lower, upper = current_upper, tol = .Machine$double.eps)$root

      success_uniroot <- TRUE
      break
    }, error = function(e) {
      # expand the search space if the saddlepoint is not found
      expansion_factor <- ifelse(i <= max_expansions/2, 2, 10)
      current_lower <<- current_lower * expansion_factor
      current_upper <<- current_upper * expansion_factor
    }
    )

    if(success_uniroot == TRUE) break
  }

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
                   gcm.default = FALSE)
  }else{
    test_stat <- sum(prod_resids)/(stats::sd(prod_resids) * sqrt(n-1))

    res <- list(test_stat = test_stat,
                p.left = stats::pnorm(test_stat, lower.tail = TRUE),
                p.right = stats::pnorm(test_stat, lower.tail = FALSE),
                p.both = 2*stats::pnorm(abs(test_stat), lower.tail = FALSE),
                gcm.default = TRUE)
  }

  return(res)
}



#####################################################################################
#' \code{dCRT_dist} is a function that returns simulated data from an appropriate
#' distribution depending on a specified  GLM family
#'
#' @param n The point where the CGF will be computed.
#' @param fitted.val A vector containing the fitted parameter values by
#' fitting a GLM to X on Z.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @return Simulated data from an appropriate distribution.
#' @keywords internal
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


#####################################################################################
#' \code{wcgf} is a function computing the cumulant generating function (CGF) of
#' distributions, multiplied by a weight function, from the GLM family
#'
#' @param s The point where the CGF will be computed.
#' @param P A vector containing the parameter values of the family of distributions.
#' @param W A vector containing the weights.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @return CGF of the weighted distribution evaluated at \code{s}.
#' @keywords internal
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


#####################################################################################
#' \code{d1_wcgf} is a function computing the derivative of the cumulant generating
#' function (CGF) of distributions, multiplied by a weight function, from GLM family
#'
#' @param s The point where the CGF will be computed.
#' @param P A vector containing the parameter values of the family of distributions.
#' @param W A vector containing the weights.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @return CGF of the weighted distribution evaluated at \code{s}.
#' @keywords internal
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


####################################################################################
#' \code{d2_wcgf} is a function computing the hessian of the cumulant generating
#' function (CGF) of distributions, multiplied by a weight function, from GLM family
#'
#' @param s The point where the CGF will be computed.
#' @param P A vector containing the parameter values of the family of distributions.
#' @param W A vector containing the weights.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @return CGF of the weighted distribution evaluated at \code{s}.
#' @keywords internal
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


#######################################################################################################
#' \code{nb_precomp} is a function computing the dispersion parameter in negative binomial regression
#'
#' @param data A list containing the response Y and covariate Z
#' @return a list containing the Poisson model fitted values and estimate for dispersion
#' @keywords internal
nb_precomp <- function(data){

  Y <- data$Y; Z <- data$Z

  pois_fit <- stats::glm.fit(y = Y, x = Z, family = stats::poisson())

  theta_hat <- estimate_theta(
    y = Y,
    mu = pois_fit$fitted.values,
    dfr = pois_fit$df.residual,
    limit = 50,
    eps = (.Machine$double.eps)^(1/4)
  )[[1]]

  return(list(fitted_values = pois_fit$fitted.values, theta_hat = theta_hat))
}




#####################################################################################
#' Fit models to data
#'
#' \code{fit_models} is a function carrying out the saddlepoint approximation to the
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
#' @param fit_vals_own Fitted values of X_on_Z and Y_on_Z based on user's own method
#' @return Simulated data from an appropriate distribution.
#' @keywords internal
fit_models <- function(data,
                       X_on_Z_fam, Y_on_Z_fam,
                       fitting_method = 'glm',
                       fit_vals_own = NULL){

   # extract (X,Y,Z) from inputted data
   X <- data$X; Y <- data$Y; Z <- data$Z
   additional_info <- list()

   if(fitting_method == 'glm'){
      # fit X on Z regression when fitting method is glm
      X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = X_on_Z_fam))
      X_on_Z_fit_vals <- X_on_Z_fit$fitted.values

      # fit Y on Z regression when fitting method is glm
      if(Y_on_Z_fam == "negative.binomial"){
         aux_info_Y_on_Z <- nb_precomp(list(Y = Y, Z = Z))

         Y_on_Z_fit <- stats::glm(Y ~ Z,
                                  family = MASS::negative.binomial(aux_info_Y_on_Z$theta_hat),
                                  mustart = aux_info_Y_on_Z$fitted_values) |> suppressWarnings()
         Y_on_Z_fit_vals <- Y_on_Z_fit$fitted.values

         additional_info$NB.disp.param <- aux_info_Y_on_Z$theta_hat
      } else{
         Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = Y_on_Z_fam))
         Y_on_Z_fit_vals <- Y_on_Z_fit$fitted.values

         additional_info$NB.disp.param <- NA
      }
   } else if(fitting_method == 'rf' || fitting_method == "prob_forest"){
      # fit X on Z regression when fitting method is random forest
      if(X_on_Z_fam == "binomial"){
         p.forest.X <- grf::probability_forest(X = as.matrix(Z), Y = as.factor(X))
         p.hat.X <- stats::predict(p.forest.X, as.matrix(Z), estimate.variance = F)

         X_on_Z_fit_vals <- p.hat.X$predictions[ ,"1"]
      }

      # fit Y on Z regression when fitting method is random forest
      if(Y_on_Z_fam == "binomial"){
         p.forest.Y <- grf::probability_forest(X = as.matrix(Z), Y = as.factor(Y))
         p.hat.Y <- stats::predict(p.forest.Y, as.matrix(Z), estimate.variance = F)

         Y_on_Z_fit_vals <- p.hat.Y$predictions[ ,"1"]
      }

      additional_info$NB.disp.param <- NA
   } else if(fitting_method == 'own') {
      # Validate that fit_vals_own is provided and contains necessary components
      if(is.list(fit_vals_own) ||
          !all(c("X_on_Z_fit_vals", "Y_on_Z_fit_vals") %in% names(fit_vals_own))) {
         stop("fit_vals_own must be a list containing 'X_on_Z_fit_vals' and 'Y_on_Z_fit_vals' when using fitting_method = 'own'")
      }

      # Extract pre-computed values from fit_vals_own
      X_on_Z_fit_vals <- fit_vals_own$X_on_Z_fit_vals
      Y_on_Z_fit_vals <- fit_vals_own$Y_on_Z_fit_vals

      additional_info$NB.disp.param <- NA
   }

   return(list(X_on_Z_fit_vals = X_on_Z_fit_vals,
               Y_on_Z_fit_vals = Y_on_Z_fit_vals,
               additional_info = additional_info))
}








# spacrt - method_helpers.R
