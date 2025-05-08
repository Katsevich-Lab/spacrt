library(testthat)

set.seed(1903)

n <- 100; p <- 6

data <- list(X = rbinom(n = n, size = 1, prob = 0.4),
             Y = rpois(n = n, lambda = 2),
             Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))

X <- data$X; Y <- data$Y; Z <- data$Z

X_on_Z_fam <- "binomial"; Y_on_Z_fam <- "poisson"
fitting_X_on_Z <- 'rf'; fitting_Y_on_Z <- 'glm'

results_GCM <- GCM(X, Y, Z,
                   family = list(XZ = X_on_Z_fam, YZ = Y_on_Z_fam),
                   method = list(XZ = fitting_X_on_Z, YZ = fitting_Y_on_Z))

results_dCRT <- dCRT(X, Y, Z,
                     family = list(XZ = X_on_Z_fam, YZ = Y_on_Z_fam),
                     method = list(XZ = fitting_X_on_Z, YZ = fitting_Y_on_Z),
                     B = 10000)

results_spaCRT <- spaCRT(X, Y, Z,
                         family = list(XZ = X_on_Z_fam, YZ = Y_on_Z_fam),
                         method = list(XZ = fitting_X_on_Z, YZ = fitting_Y_on_Z))

results_test_stat <- c(results_GCM$test_stat, results_dCRT$test_stat, results_spaCRT$test_stat)
results_p_value <- c(results_GCM$p_value, results_dCRT$p_value, results_spaCRT$p_value)


test_that("All test statistics and p-values were accurately computed.",{
  expect_equal(object = results_test_stat |> round(5),
               expected = c(0.05634,  0.00019, -0.00844))
  expect_equal(object = results_p_value |> unname() |> round(6),
               expected = c(0.955067, 0.981902, 0.986141 ))
  expect_equal(object = results_spaCRT$spa.success,
               expected = TRUE)
})

