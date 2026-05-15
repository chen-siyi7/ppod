test_that("simulate_trial returns a valid result on S4b", {
  set.seed(1)
  scen <- ppod_scenarios$S4b
  res <- simulate_trial(scen, criterion = "joint",
                         N_trial = 18L, cohort = 3L, n_runin = 9L,
                         S_interim = 50L, M_interim = 20L,
                         S_final = 100L, M_final = 30L)
  expect_true(is.list(res))
  expect_true(all(c("selected", "dose_history", "Y_T", "Y_E",
                    "stopped_tox", "stopped_fut") %in% names(res)))
  if (!is.na(res$selected)) {
    expect_gte(res$selected, 1L)
    expect_lte(res$selected, length(scen$mu_T))
  }
  expect_length(res$Y_T, length(res$Y_E))
  expect_length(res$Y_T, length(res$dose_history))
})

test_that("simulate_trial respects no-skipping rule", {
  set.seed(2)
  scen <- ppod_scenarios$S4b
  res <- simulate_trial(scen, criterion = "joint",
                         N_trial = 27L, cohort = 3L, n_runin = 9L,
                         S_interim = 50L, M_interim = 20L,
                         S_final = 100L, M_final = 30L)
  hist <- res$dose_history
  # Successive dose changes by at most +1 in each cohort step.
  cohort_starts <- seq(1L, length(hist), by = 3L)
  doses_at_cohort <- hist[cohort_starts]
  jumps <- diff(doses_at_cohort)
  expect_true(all(jumps <= 1L))
})

test_that("boin12_trial completes and returns expected fields", {
  set.seed(3)
  scen <- ppod_scenarios$S4b
  res <- boin12_trial(scen, N_trial = 18L, cohort = 3L)
  expect_true(all(c("selected", "dose_history", "n_dose", "nT", "nE",
                    "stopped_tox", "stopped_fut") %in% names(res)))
  expect_length(res$n_dose, length(scen$mu_T))
  expect_length(res$nT, length(scen$mu_T))
  expect_length(res$nE, length(scen$mu_T))
  expect_equal(sum(res$n_dose), length(res$dose_history))
})

test_that("run_replicates returns a list of the expected length", {
  scen <- ppod_scenarios$S4b
  cell <- run_replicates(scen, design = "joint", n_reps = 3L, seed_base = 100L,
                          N_trial = 18L, cohort = 3L, n_runin = 9L,
                          S_interim = 50L, M_interim = 20L,
                          S_final = 100L, M_final = 30L)
  expect_length(cell, 3L)
  expect_true(all(vapply(cell, is.list, logical(1L))))
})

test_that("operating_characteristics handles a small cell", {
  set.seed(4)
  scen <- ppod_scenarios$S4b
  cell <- run_replicates(scen, design = "joint", n_reps = 5L, seed_base = 200L,
                          N_trial = 18L, cohort = 3L, n_runin = 9L,
                          S_interim = 50L, M_interim = 20L,
                          S_final = 100L, M_final = 30L)
  oc <- operating_characteristics(cell, scen)
  expect_true(is.list(oc))
  expect_true(all(c("p_correct", "se_correct", "p_near", "p_no_select",
                    "mean_loss", "selection_table", "alloc",
                    "mean_tox_burden") %in% names(oc)))
  expect_true(oc$p_correct >= 0 && oc$p_correct <= 1)
})
