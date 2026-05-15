test_that("ppod_scenarios contains 12 named scenarios", {
  expected_names <- c("S1", "S2", "S3", "S4a", "S4b", "S4c", "S4d",
                      "S5", "S6", "S7a", "S7b", "S7c")
  expect_equal(names(ppod_scenarios), expected_names)
})

test_that("each scenario has the required fields", {
  required <- c("mu_T", "mu_E", "rho", "sigma_T", "sigma_E", "util",
                "pi_T", "pi_E", "true_OBD")
  for (name in names(ppod_scenarios)) {
    scen <- ppod_scenarios[[name]]
    for (field in required) {
      expect_true(field %in% names(scen),
                  info = paste(name, "missing", field))
    }
    expect_equal(length(scen$mu_T), length(scen$mu_E),
                 info = paste(name, "mu_T/mu_E length mismatch"))
    expect_true(length(scen$rho) == 1L ||
                length(scen$rho) == length(scen$mu_T),
                info = paste(name, "rho length mismatch"))
  }
})

test_that("make_scenario auto-computes true_OBD when not supplied", {
  s <- make_scenario(
    mu_T = c(-2, -1, 0, 1, 2),
    mu_E = c(-1, 0, 1, 0, -1),     # peak at d_3
    rho  = 0
  )
  expect_equal(s$true_OBD, 3L)
})

test_that("S4d rho_j is dose-varying as specified", {
  expect_equal(ppod_scenarios$S4d$rho,
               c(-0.85, -0.85, 0.95, 0.92, 0.9))
})
