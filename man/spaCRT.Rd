% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{spaCRT}
\alias{spaCRT}
\title{The saddlepoint approximation to the dCRT.}
\usage{
spaCRT(
  data,
  X_on_Z_fam,
  Y_on_Z_fam,
  fitting_X_on_Z = "glm",
  fitting_Y_on_Z = "glm",
  R = 5
)
}
\arguments{
\item{data}{A named list with fields \code{X} (an nx1 vector for the predictor
variable of interest), \code{Y} (an nx1 response vector), and \code{Z}
(an nxp matrix of covariates).}

\item{X_on_Z_fam}{The GLM family for the regression of X on Z
(values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).}

\item{Y_on_Z_fam}{The GLM family for the regression of Y on Z
(values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).}

\item{fitting_X_on_Z}{The fitting method for the regression X on Z.}

\item{fitting_Y_on_Z}{The fitting method for the regression Y on Z.}

\item{R}{Upper bound of search space for the saddlepoint}
}
\value{
A named list with fields \code{test_stat}, \code{left_side_p_value},
\code{right_side_p_value}, \code{both_side_p_value}, and
\code{gcm.default}.
gcm.default returns TRUE if GCM was employed due to the failure of spaCRT.
}
\description{
\code{spaCRT} is a function carrying out the saddlepoint approximation to the
dCRT based on GLMs for \code{X|Z} and \code{Y|Z}.
}
\examples{
n <- 50; p <- 4
data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
             Y = rbinom(n = n, size = 1, prob = 0.7),
             Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
X_on_Z_fam <- "binomial"
Y_on_Z_fam <- "binomial"
spaCRT(data, X_on_Z_fam, Y_on_Z_fam,
       fitting_X_on_Z = 'rf',
       fitting_Y_on_Z = 'glm')

}
