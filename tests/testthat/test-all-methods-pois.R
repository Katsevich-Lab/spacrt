library(testthat)

set.seed(2203)

n <- 100; p <- 6

data <- list(X = rbinom(n = n, size = 1, prob = 0.4),
             Y = rpois(n = n, lambda = 1),
             Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))

X_on_Z_fam <- "binomial"; Y_on_Z_fam <- "poisson"
fitting_X_on_Z <- 'glm'; fitting_Y_on_Z <- 'glm'

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
               expected = c(-0.38643, -0.15740, -0.15740))
  expect_equal(object = results_p_value |> unname() |> round(6),
               expected = c(0.699180, 0.723928, 0.727422))
  expect_equal(object = results_spaCRT$spa.success,
               expected = TRUE)
})

