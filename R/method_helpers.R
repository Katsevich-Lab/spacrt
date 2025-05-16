#' @useDynLib spacrt, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


############################################################################################
#' \code{spa_cdf} SPA to CDF of T_n = S_n / sqrt(n)
#'
#' @param X Numeric vector of length n for the predictor variable.
#' @param Y Numeric vector of length n for the response variable.
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
#'                  fam = "binomial", R = 100)
#'
#' @keywords internal
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

  success_uniroot <- FALSE

  tryCatch({
    # try to solve the saddlepoint equation
    s.hat <- stats::uniroot(
      f = function(s) d1_wcgf(s, P = P, W = W, fam) - sqrt(n)*t,
      lower = -abs(R), upper = abs(R),
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
#' @param n Number of samples
#' @param fitted.val A vector containing the fitted parameter values by
#' fitting a GLM to X on Z.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#'
#' @return Simulated data from an appropriate distribution.
#'
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
#'
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


############################################################################################
#' \code{d1_wcgf} is a function computing the derivative of the weighted cumulant
#' generating function (WCGF) of distributions from GLM family
#'
#' @inheritParams wcgf
#'
#' @return First derivative of CGF of the weighted distribution evaluated at \code{s}.
#'
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


############################################################################################
#' \code{d2_wcgf} is a function computing the hessian of the weighted cumulant
#' generating function (WCGF) of distributions from GLM family
#'
#' @inheritParams wcgf
#'
#' @return Second derivative of CGF of the weighted distribution evaluated at \code{s}.
#'
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


############################################################################################
#' \code{nb_precomp} is a function computing the dispersion parameter in negative
#' binomial regression
#'
#' @param V A list containing the response
#' @param Z A list containing the covariate
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{fitted_values}{The fitted values from the Poisson regression of \code{Y} on \code{Z}.}
#'   \item{theta_hat}{The estimated dispersion parameter (theta) for the negative binomial model, computed via maximum likelihood or method of moments.}
#' }
#'
#' @keywords internal
nb_precomp <- function(V,Z){

  # Fit Poisson GLM: Y ~ Z with intercept
  pois_fit <- stats::glm.fit(y = V, x = cbind(1, Z), family = stats::poisson())

  # Estimate NB dispersion parameter using C++ function estimate_theta
  theta_hat <- estimate_theta(
    y = V,
    mu = pois_fit$fitted.values,
    dfr = pois_fit$df.residual,
    limit = 200,
    eps = (.Machine$double.eps)^(1/2)
  )[[1]]

  # if truncate theta at 5e-3 if it is less than that
  theta_hat <- max(5e-3, theta_hat)

  return(list(fitted_values = pois_fit$fitted.values,
              theta_hat = theta_hat))
}


############################################################################################
#' \code{fit_models} is a function carrying out the regressions `X` on `Z` and `Y` on `Z`.
#'
#' @inheritParams GCM
#'
#' @return A named list of fitted values of X|Z and Y|Z.
#'
#' @keywords internal
fit_models <- function(X, Y, Z,
                       family,
                       method = list(XZ = 'glm', YZ = 'glm'),
                       fitted = list(XZ = NULL, YZ = NULL)){

   # fit X on Z regression
   X_on_Z_fit_vals <- fit_single_model(V = X, Z = Z,
                                       V_on_Z_fam = family$XZ,
                                       fitting_V_on_Z = method$XZ,
                                       fit_vals_V_on_Z_own = fitted$XZ)

   # fit Y on Z regression
   Y_on_Z_fit_vals <- fit_single_model(V = Y, Z = Z,
                                       V_on_Z_fam = family$YZ,
                                       fitting_V_on_Z = method$YZ,
                                       fit_vals_V_on_Z_own = fitted$YZ)

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
#' (values can be \code{glm} (default) or \code{random_forest}).
#' @param fit_vals_V_on_Z_own Vector of fitted values for V on Z in case the user's custom method.
#' Works only if fitting_V_on_Z = 'own'.
#'
#' @return A vector of fitted values of V|Z.
#'
#' @keywords internal
fit_single_model <- function(V, Z,
                             V_on_Z_fam,
                             fitting_V_on_Z,
                             fit_vals_V_on_Z_own){

  if(!is.null(fitting_V_on_Z)){
    if(fitting_V_on_Z == 'glm'){
      # fit V on Z regression when fitting method is glm
      if(V_on_Z_fam == "negative.binomial"){
        # V_on_Z_fam == "negative.binomial"
        aux_info_V_on_Z <- nb_precomp(V = V, Z = Z)

        V_on_Z_fit <- stats::glm(V ~ Z,
                                 family = MASS::negative.binomial(aux_info_V_on_Z$theta_hat),
                                 mustart = aux_info_V_on_Z$fitted_values) |> suppressWarnings()

        V_on_Z_fit_vals <- V_on_Z_fit$fitted.values
      } else{
        # V_on_Z_fam == any other glm family
        V_on_Z_fit <- suppressWarnings(stats::glm(V ~ Z, family = V_on_Z_fam))
        V_on_Z_fit_vals <- V_on_Z_fit$fitted.values
      }
    } else if(fitting_V_on_Z == 'random_forest'){
      # fit V on Z regression when fitting method is random forest
      Z <- as.data.frame(Z)
      colnames(Z) <- paste0("V", seq_len(ncol(Z)))

      discrete_fam <- c('binomial','poisson','negative.binomial')

      if(V_on_Z_fam %in% discrete_fam){
        # V_on_Z_fam is discrete
        rf_fit <- ranger::ranger(y = as.factor(V), x = Z, probability = TRUE)
        pred_probs <- stats::predict(rf_fit, data = Z)$predictions

        # Extract probability for class "1" if present; otherwise use first level
        target_class <- if("1" %in% colnames(pred_probs)) "1" else colnames(pred_probs)[1]
        V_on_Z_fit_vals <- pred_probs[, target_class]
      } else{
        # V_on_Z_fam is continuous
        rf_fit <- ranger::ranger(y = V, x = Z)
        V_on_Z_fit_vals <- stats::predict(rf_fit, data = Z)$predictions
      }
    }
  } else{
    # Validate that fit_vals_V_on_Z_own is provided and contains necessary components
    if(!is.numeric(fit_vals_V_on_Z_own)){
      stop("fit_vals_V_on_Z_own must be a vector containing 'V_on_Z_fit_vals' when using fitting_V_on_Z = 'own'")
    }

    # Extract pre-computed values from fit_vals_V_on_Z_own
    V_on_Z_fit_vals <- fit_vals_V_on_Z_own
  }

  return(V_on_Z_fit_vals)
}

############################################################################################
#' Validate inputs to the GCM function
#'
#' \code{check_inputs_main} is a utility function that checks whether the arguments provided to
#' \code{GCM()}, \code{dCRT()} and \code{spaCRT()} satisfy the required structure and types.
#'
#' @param X Numeric vector of length \eqn{n}, representing the predictor variable.
#' @param Y Numeric vector of length \eqn{n}, representing the response variable.
#' @param Z Numeric matrix with \eqn{n} rows and \eqn{p} columns, representing covariates.
#' @param family Named list with elements \code{XZ} and \code{YZ}, each a character string
#'   specifying the model family for \eqn{X \mid Z} and \eqn{Y \mid Z}, respectively.
#' @param method Named list with elements \code{XZ} and \code{YZ}, each a character string
#'   indicating the model-fitting method to use, e.g., \code{"glm"}, \code{"random_forest"}.
#' @param fitted Named list with elements \code{XZ} and \code{YZ}, each either \code{NULL}
#'   or a numeric vector of length \eqn{n} representing user-supplied fitted values.
#' @param alternative Character string indicating the alternative hypothesis for the test.
#'   Must be one of \code{"two.sided"}, \code{"greater"}, or \code{"less"}.
#'
#' @return None. The function throws an error if any check fails.
#'
#' @keywords internal
check_inputs_main <- function(X, Y, Z,
                              family, method,
                              fitted,
                              alternative,
                              func = 'GCM', B = NULL) {
  n <- length(X)

  # Check X
  if(!is.numeric(X) || !is.vector(X)) stop("`X` must be a numeric vector.")

  # Check Y
  if(!is.numeric(Y) || !is.vector(Y)) stop("`Y` must be a numeric vector.")
  if(length(Y) != n) stop("`X` and `Y` must be of the same length.")

  # Check Z
  if(!is.matrix(Z) || !is.numeric(Z)) stop("`Z` must be a numeric matrix.")
  if(nrow(Z) != n) stop("`Z` must have the same number of rows as `X` and `Y`.")

  # supported families and fitting methods
  supported_families <- c("binomial", "poisson","negative.binomial")
  supported_methods  <- c("glm", "random_forest")

  # Check fitted
  if (!is.list(fitted) || !setequal(names(fitted), c("XZ", "YZ"))) {
    stop("`fitted` must be a named list with exactly two elements: 'XZ' and 'YZ'.")
  }
  for (name in c("XZ", "YZ")) {
    val <- fitted[[name]]
    if (!is.null(val)) {
      if (!is.numeric(val) || !is.vector(val) || length(val) != n)
        stop(sprintf("`fitted$%s` must be NULL or a numeric vector of length %d.", name, n))
    }
  }

  # Check family
  if(!is.list(family)) stop("`family` must be a named list.")

  for (name in c("XZ", "YZ")) {
    if (is.null(fitted[[name]])) {
      # Must be provided and valid
      if (!(name %in% names(family)))
        stop(sprintf("`family` must include a character entry for '%s' if `fitted$%s` is NULL.", name, name))
      family_val <- family[[name]]
      if (!is.character(family_val) || length(family_val) != 1)
        stop(sprintf("`family$%s` must be a character string.", name))
      if (!(family_val %in% supported_families))
        stop(sprintf("`family$%s` must be one of: %s.", name, paste(supported_families, collapse = ", ")))
    }
  }

  # Check method
  if(!is.list(method)) stop("`method` must be a named list.")

  for (name in c("XZ", "YZ")) {
    if (is.null(fitted[[name]])) {
      # Must be provided and valid
      if (!(name %in% names(method)))
        stop(sprintf("`method` must include a character entry for '%s' if `fitted$%s` is NULL.", name, name))
      method_val <- method[[name]]
      if (!is.character(method_val) || length(method_val) != 1)
        stop(sprintf("`method$%s` must be a character string.", name))
      if (!(method_val %in% supported_methods))
        stop(sprintf("`method$%s` must be one of: %s.", name, paste(supported_methods, collapse = ", ")))
    }
  }

  # Check alternative
  if(!alternative %in% c("two.sided", "greater", "less"))
    stop("`alternative` must be one of 'two.sided', 'greater', or 'less'.")

  # Special check for dCRT
  if (func == 'dCRT') {
    if (!is.numeric(B) || length(B) != 1 || B <= 0 || B != as.integer(B)) {
      stop("`B` must be a positive integer.")
    }
  }

  return(invisible(NULL))
}


# spacrt - method_helpers.R
