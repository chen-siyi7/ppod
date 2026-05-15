test_that("make_utility returns a function with the documented form", {
  u <- make_utility(alpha = 1, beta = 1, lam = 2)
  expect_type(u, "closure")
  # At x_T = 0, x_E = 1, utility = 1
  expect_equal(u(x_T = 0, x_E = 1), 1)
  # At x_T = 0, x_E = 0, utility = 0
  expect_equal(u(x_T = 0, x_E = 0), 0)
  # At x_T = 1, x_E = 1, utility = exp(-lam)
  expect_equal(u(x_T = 1, x_E = 1), exp(-2))
})

test_that("make_utility validates parameter ranges", {
  expect_error(make_utility(alpha = 0.5, beta = 1, lam = 2))
  expect_error(make_utility(alpha = 1, beta = 0.5, lam = 2))
  expect_error(make_utility(alpha = 1, beta = 1, lam = -1))
})

test_that("utility_substitution matches make_utility", {
  u <- make_utility(1, 2, 3)
  expect_equal(
    u(x_T = 0.2, x_E = 0.7),
    utility_substitution(x_T = 0.2, x_E = 0.7, alpha = 1, beta = 2, lam = 3)
  )
})

test_that("utility is vectorized", {
  u <- make_utility(1, 1, 2)
  result <- u(x_T = c(0.1, 0.3, 0.5), x_E = c(0.5, 0.7, 0.9))
  expect_length(result, 3L)
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

test_that("utility has negative mixed partial (substitution condition)", {
  # U_TE < 0: numerically verify on a grid
  u <- make_utility(1, 1, 2)
  eps <- 1e-4
  xs <- seq(0.1, 0.9, by = 0.2)
  for (xT in xs) {
    for (xE in xs) {
      u_pp <- u(xT + eps, xE + eps)
      u_pm <- u(xT + eps, xE - eps)
      u_mp <- u(xT - eps, xE + eps)
      u_mm <- u(xT - eps, xE - eps)
      cross <- (u_pp - u_pm - u_mp + u_mm) / (4 * eps^2)
      expect_lt(cross, 0)
    }
  }
})
