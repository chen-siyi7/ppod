test_that("true_expected_utility matches S4d truth-table values from the manuscript", {
  scen <- ppod_scenarios$S4d
  eu_joint_truth <- true_expected_utility(scen, n_mc = 200000L)
  # Manuscript Table 2: Psi_j(rho_j) approximately
  #   d_1: 0.1899, d_2: 0.3461, d_3: 0.3267, d_4: 0.2290, d_5: 0.1252
  # MCSE per dose ~ 0.0002 at 2e6 draws; with 2e5 draws, MCSE ~ 0.0006
  # Allow 0.005 tolerance to absorb MC noise
  expected <- c(0.1899, 0.3461, 0.3267, 0.2290, 0.1252)
  expect_equal(eu_joint_truth, expected, tolerance = 0.005)
})

test_that("marginal vs joint true OBD reverses in S4d", {
  scen <- ppod_scenarios$S4d
  eu_j <- true_expected_utility(scen, n_mc = 200000L)
  eu_m <- true_expected_utility(scen, n_mc = 200000L, marginal = TRUE)
  # Joint OBD is d_2; marginal OBD is d_3
  expect_equal(which.max(eu_j), 2L)
  expect_equal(which.max(eu_m), 3L)
})

test_that("S7a joint utility favors d_3 by a small margin", {
  scen <- ppod_scenarios$S7a
  eu_j <- true_expected_utility(scen, n_mc = 500000L)
  # Manuscript: Psi_2(rho) = 0.3946, Psi_3(rho) = 0.3963
  # gap of +0.0017 in favor of d_3
  expect_gt(eu_j[3] - eu_j[2], 0)
  expect_lt(eu_j[3] - eu_j[2], 0.005)
})

test_that("S1 through S6 true OBD matches the documented true_OBD", {
  for (key in c("S1", "S2", "S3", "S4a", "S4b", "S4c", "S5", "S6")) {
    scen <- ppod_scenarios[[key]]
    eu <- true_expected_utility(scen, n_mc = 100000L)
    expect_equal(which.max(eu), scen$true_OBD,
                 info = paste("scenario", key))
  }
})
