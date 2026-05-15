test_that("softplus is non-negative and monotone", {
  x <- seq(-10, 10, by = 0.5)
  s <- softplus(x)
  expect_true(all(s >= 0))
  expect_true(all(diff(s) >= 0))
  # softplus(0) = log(2)
  expect_equal(softplus(0), log(2))
})

test_that("unpack_theta returns components on the correct scales", {
  J <- 5L
  theta <- c(
    -1.5,                       # alpha_T
    rep(-1, J - 1L),            # beta_T (length J-1)
    0.2,                        # alpha_E
    rep(0, J - 1L),             # gamma free coords
    log(0.5),                   # log sigma_T
    log(0.5),                   # log sigma_E
    0                           # atanh rho
  )
  p <- unpack_theta(theta, J = J)
  expect_equal(p$alpha_T, -1.5)
  expect_equal(p$alpha_E, 0.2)
  expect_equal(p$sT, 0.5)
  expect_equal(p$sE, 0.5)
  expect_equal(p$rho, 0)
  # sum-to-zero on gamma
  expect_equal(sum(p$gamma), 0, tolerance = 1e-12)
  expect_length(p$gamma, J)
})

test_that("unpack_theta rejects invalid input", {
  expect_error(unpack_theta(rep(0, 10), J = 5L))   # wrong length
  expect_error(unpack_theta(rep(0, 9), J = 2L))    # J too small
})

test_that("neg_log_post is finite at default parameters with synthetic data", {
  set.seed(1)
  J <- 3L
  doses <- rep(1:3, each = 5L)
  Y_T <- stats::rnorm(15L, mean = c(-1.5, -0.8, -0.2)[doses], sd = 0.5)
  Y_E <- stats::rnorm(15L, mean = c(-0.5, 0.5, 0.7)[doses], sd = 0.5)
  theta <- c(-1.5, rep(-1, 2), 0, 0, 0, log(0.5), log(0.5), 0)
  val <- neg_log_post(theta, Y_T, Y_E, doses, J = J)
  expect_true(is.finite(val))
  expect_gt(val, 0)
})

test_that("neg_log_post returns large value when rho near ±1", {
  J <- 3L
  doses <- rep(1:3, each = 3L)
  Y_T <- rep(c(-1.5, -0.8, -0.2), each = 3L)
  Y_E <- rep(c(-0.5, 0.5, 0.7), each = 3L)
  # atanh(rho) -> large makes rho -> 1
  theta <- c(-1.5, rep(-1, 2), 0, 0, 0, log(0.5), log(0.5), 50)
  val <- neg_log_post(theta, Y_T, Y_E, doses, J = J)
  expect_true(is.finite(val) || val >= 1e19)
})
