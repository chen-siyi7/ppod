test_that("fit_laplace converges on well-conditioned synthetic data", {
  set.seed(20260515)
  J <- 5L
  n_per <- 20L
  doses <- rep(seq_len(J), each = n_per)
  true_mu_T <- c(-2.2, -1.5, -0.8, -0.2, 0.4)
  true_mu_E <- c(-1.4, -0.2, 0.6, 1.4, 2.2)
  Y_T <- stats::rnorm(J * n_per, mean = true_mu_T[doses], sd = 0.5)
  Y_E <- stats::rnorm(J * n_per, mean = true_mu_E[doses], sd = 0.5)

  fit <- fit_laplace(Y_T, Y_E, doses, J = J)
  expect_true(fit$ok)
  expect_length(fit$theta, 2L * J + 3L)
  expect_equal(dim(fit$cov), c(2L * J + 3L, 2L * J + 3L))
  # Covariance is positive-definite
  ev <- eigen(fit$cov, only.values = TRUE)$values
  expect_true(all(ev > 0))
})

test_that("fit_laplace recovers dose-level means approximately at large n", {
  set.seed(1)
  J <- 5L
  n_per <- 60L
  doses <- rep(seq_len(J), each = n_per)
  true_mu_T <- c(-2.2, -1.5, -0.8, -0.2, 0.4)
  true_mu_E <- c(-1.4, -0.2, 0.6, 1.4, 2.2)
  Y_T <- stats::rnorm(J * n_per, mean = true_mu_T[doses], sd = 0.5)
  Y_E <- stats::rnorm(J * n_per, mean = true_mu_E[doses], sd = 0.5)
  fit <- fit_laplace(Y_T, Y_E, doses, J = J)
  p <- unpack_theta(fit$theta, J = J)
  # Reconstruct mu_T from softplus increments
  mu_T_hat <- p$alpha_T + c(0, cumsum(softplus(p$beta_T)))
  mu_E_hat <- p$alpha_E + p$gamma
  expect_equal(mu_T_hat, true_mu_T, tolerance = 0.15)
  expect_equal(mu_E_hat, true_mu_E, tolerance = 0.25)
})

test_that("sample_posterior draws have the requested shape and approximate covariance", {
  set.seed(2)
  J <- 5L
  n_per <- 30L
  doses <- rep(seq_len(J), each = n_per)
  Y_T <- stats::rnorm(J * n_per, mean = c(-2, -1, 0, 1, 2)[doses], sd = 0.5)
  Y_E <- stats::rnorm(J * n_per, mean = c(-1, 0, 1, 0.5, 0)[doses], sd = 0.5)
  fit <- fit_laplace(Y_T, Y_E, doses, J = J)
  S <- 2000L
  draws <- sample_posterior(fit$theta, fit$cov, S = S)
  expect_equal(dim(draws), c(S, 2L * J + 3L))
  # Sample mean ~ theta_hat
  expect_equal(colMeans(draws), fit$theta, tolerance = 0.15)
  # Sample covariance ~ fit$cov (loose tolerance)
  emp_cov <- cov(draws)
  expect_equal(diag(emp_cov), diag(fit$cov), tolerance = 0.2)
})
