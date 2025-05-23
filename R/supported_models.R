#' Supported Families and Model-fitting Methods
#'
#' This page documents the families and fitting methods supported in the \code{X | Z} and
#' \code{Y | Z} models for functions \code{\link{GCM}}, \code{\link{dCRT}}, and \code{\link{spaCRT}}.
#'
#' \strong{Supported Families for \code{X | Z}:}
#' \itemize{
#'   \item \code{"binomial"}: Models binary-valued (0/1) \code{X} using logistic regression
#'     or classification forests, depending on the fitting method specified.
#'   \item \code{"poisson"}: Models count-valued \code{X} under a Poisson distribution,
#'     with the canonical log link. Appropriate for non-negative integer outcomes.
#' }
#'
#' Other families may be added in future versions.
#'
#' \strong{Supported Families for \code{Y | Z}:}
#' \itemize{
#'   \item \code{"binomial"}: Models binary-valued (0/1) \code{Y} using logistic regression
#'     or classification forests, depending on the fitting method specified.
#'   \item \code{"poisson"}: Models count-valued \code{Y} under a Poisson distribution,
#'     with the canonical log link. Appropriate for non-negative integer outcomes.
#'   \item \code{"negative.binomial"}: Used when \code{Y} exhibits overdispersed
#'     count behavior that cannot be captured by a Poisson model. The dispersion
#'     parameter (theta) is estimated internally using a vendored `C++` implementation
#'     based on code from the \pkg{sceptre} package. The procedure first attempts
#'     maximum likelihood estimation via Newton's method, with a fallback to a
#'     method-of-moments estimator if convergence fails.
#' }
#'
#' Other families may be added in future versions.
#'
#' \strong{Supported Model-fitting Methods:}
#' \itemize{
#'   \item \code{"glm"}: Fits a generalized linear model using \code{stats::glm()},
#'     with the model family determined by the specified family string
#'     (e.g., \code{"binomial"}, \code{"poisson"}, etc.). An intercept is automatically
#'     included in the model; users should provide \code{Z} without an intercept column.
#'   \item \code{"random_forest"}: Fits a random forest using \code{ranger::ranger()} with
#'     classification or regression mode depending on the family.
#' }
#'
#' Other methods may be added in future versions.
#'
#' @name supported_models
#' @aliases supported_models
#' @title Supported Model Families and Model-fitting Methods for Inference
#' @keywords documentation
#' @seealso \code{\link{GCM}}, \code{\link{dCRT}}, \code{\link{spaCRT}}
NULL
