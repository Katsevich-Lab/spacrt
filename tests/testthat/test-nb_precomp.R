# generated using ```usethis::use_test("nb_precomp")```

test_that("nb_precomp returns fitted values and theta_hat", {
  # test_env <- new.env()
  # test_env$nb_precomp <- get("nb_precomp", envir = asNamespace("spacrt"))

  set.seed(1)

  n <- 500; p <- 2
  gamma_0 <- -3; gamma_1 <- 1
  beta_0 <- -3; beta_1 <- 1
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

  out <- spacrt:::nb_precomp(data)

  expect_type(out, "list")
  expect_named(out, c("fitted_values", "theta_hat"))
  expect_equal(length(out$fitted_values), n)
  expect_true(is.numeric(out$theta_hat))
  expect_true(out$theta_hat > 0)
  expect_equal(out$theta_hat, 95.33628)
})
