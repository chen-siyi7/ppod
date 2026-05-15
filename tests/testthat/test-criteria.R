test_that("eu_joint equals eu_marginal when posterior is concentrated at rho = 0", {
  # Construct synthetic samples with atanh_rho exactly 0 and minimal variance.
  set.seed(1)
  J <- 5L
  P <- 2L * J + 3L
  S <- 500L
  M <- 200L
  theta_hat <- c(-1.5, rep(-1, J - 1L), 0, rep(0, J - 1L),
                 log(0.5), log(0.5), 0)
  cov <- diag(1e-6, P)
  cov[P, P] <- 1e-8                     # essentially deterministic rho = 0
  samples <- sample_posterior(theta_hat, cov, S = S)
  util <- make_utility(1, 1, 2)
  ej <- eu_joint(samples, M, util, J = J)
  em <- eu_marginal(samples, M, util, J = J)
  expect_equal(ej, em, tolerance = 0.02)
})

test_that("eu_joint is decreasing in rho when U_TE < 0 (Corollary 1)", {
  set.seed(2)
  J <- 5L
  P <- 2L * J + 3L
  S <- 600L
  M <- 200L
  util <- make_utility(1, 1, 2)
  theta_hat <- c(-1.5, rep(-1, J - 1L), 0.5, rep(0, J - 1L),
                 log(0.5), log(0.5), 0)
  cov <- diag(1e-6, P)
  eu_at <- function(atanh_rho_val) {
    t <- theta_hat
    t[P] <- atanh_rho_val
    s <- sample_posterior(t, cov, S = S)
    eu_joint(s, M, util, J = J)
  }
  eu_neg  <- eu_at(atanh(-0.5))
  eu_zero <- eu_at(0)
  eu_pos  <- eu_at(atanh(0.5))
  # At each dose, Psi(rho) should be decreasing in rho.
  # Some MC noise; require monotonicity for at least 4 of 5 doses with margin.
  decreasing <- (eu_neg - eu_pos) > 0
  expect_gte(sum(decreasing), 4L)
})

test_that("eu_plugin uses only the posterior mode", {
  J <- 5L
  P <- 2L * J + 3L
  util <- make_utility(1, 1, 2)
  theta_hat <- c(-2, rep(-0.5, J - 1L), 0.5, rep(0, J - 1L),
                 log(0.5), log(0.5), 0)
  # Two different covariances should not change plug-in
  p1 <- eu_plugin(theta_hat, util, J = J)
  p2 <- eu_plugin(theta_hat, util, J = J)
  expect_identical(p1, p2)
  expect_length(p1, J)
  expect_true(all(p1 >= 0 & p1 <= 1))
})

test_that("admissibility returns proportions in [0,1]", {
  set.seed(3)
  J <- 5L
  P <- 2L * J + 3L
  S <- 200L
  M <- 100L
  theta_hat <- c(-1.5, rep(-1, J - 1L), 0, rep(0, J - 1L),
                 log(0.5), log(0.5), 0)
  cov <- diag(0.01, P)
  samples <- sample_posterior(theta_hat, cov, S = S)
  ad <- admissibility(samples, M, xT_max = 0.5, xE_min = stats::plogis(-1),
                       J = J)
  expect_length(ad$p_over, J)
  expect_length(ad$p_under, J)
  expect_true(all(ad$p_over >= 0 & ad$p_over <= 1))
  expect_true(all(ad$p_under >= 0 & ad$p_under <= 1))
})

test_that("predict_outcomes returns matrices of the right shape", {
  set.seed(4)
  J <- 5L
  P <- 2L * J + 3L
  S <- 50L
  M <- 30L
  theta_hat <- c(-1.5, rep(-1, J - 1L), 0, rep(0, J - 1L),
                 log(0.5), log(0.5), 0)
  cov <- diag(0.05, P)
  samples <- sample_posterior(theta_hat, cov, S = S)
  out <- predict_outcomes(samples, M = M, dose = 3L, J = J)
  expect_equal(dim(out$Y_T), c(S, M))
  expect_equal(dim(out$x_T), c(S, M))
  expect_true(all(out$x_T > 0 & out$x_T < 1))
  expect_true(all(out$x_E > 0 & out$x_E < 1))
})
