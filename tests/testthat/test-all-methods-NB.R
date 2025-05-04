library(testthat)

set.seed(2203)

n <- 500; p <- 4
gamma_0 <- -3
gamma_1 <- 1
beta_0 <- -3
beta_1 <- 1
rho <- 1
theta <- 1

# define data-generating model based on the Gaussian linear model
generate_data_nb <- function(n, gamma_0, gamma_1,
                             beta_0, beta_1, rho, theta){

  expit <- function(theta)(exp(theta)/(1 + exp(theta)))

  Z <- as.matrix(rnorm(n = n, mean = 0, sd = 1))
  X <- stats::rbinom(n = n, size = 1, prob = expit(gamma_0 + gamma_1*Z))
  Y <- MASS::rnegbin(n = n, mu = exp(beta_0 + beta_1*Z + rho*X), theta = theta)

  return(list(X = X, Y = Y, Z = Z))
}

data <- generate_data_nb(n, gamma_0, gamma_1,
                         beta_0, beta_1, rho, theta)

X_on_Z_fam <- "binomial"; Y_on_Z_fam <- "negative.binomial"
fitting_X_on_Z <- 'rf'; fitting_Y_on_Z <- 'glm'


results_GCM <- GCM(data, X_on_Z_fam, Y_on_Z_fam,
                   fitting_X_on_Z = fitting_X_on_Z,
                   fitting_Y_on_Z = fitting_Y_on_Z)

results_dCRT <- dCRT(data, X_on_Z_fam, Y_on_Z_fam,
                     fitting_X_on_Z = fitting_X_on_Z,
                     fitting_Y_on_Z = fitting_Y_on_Z,
                     B = 10000)

results_spaCRT <- spaCRT(data, X_on_Z_fam, Y_on_Z_fam,
                         fitting_X_on_Z = fitting_X_on_Z,
                         fitting_Y_on_Z = fitting_Y_on_Z)

results_test_stat <- c(results_GCM$test_stat, results_dCRT$test_stat, results_spaCRT$test_stat)
results_p_value <- c(results_GCM$p_value, results_dCRT$p_value, results_spaCRT$p_value)


test_that("All test statistics and p-values were accurately computed.",{
  expect_equal(object = results_test_stat |> round(5),
               expected = c(2.04561, 0.53143, 0.53400))
  expect_equal(object = results_p_value |> unname() |> round(6),
               expected = c(0.040795, 0.008799, 0.007898))
  expect_equal(object = results_spaCRT$spa.success,
               expected = TRUE)
})



