#' Run many replicates of a single (scenario, design) cell
#'
#' Runs \code{n_reps} replicates of either \code{\link{simulate_trial}}
#' (for the continuous designs \code{"joint"}, \code{"marginal"},
#' \code{"plugin"}) or \code{\link{boin12_trial}} (for
#' \code{design = "boin12"}) under a single scenario. Replicates are
#' seeded reproducibly from \code{seed_base} so that paired
#' within-replicate comparisons across designs are possible.
#'
#' @param scen Scenario list.
#' @param design One of \code{"joint"}, \code{"marginal"},
#'   \code{"plugin"}, \code{"boin12"}.
#' @param n_reps Number of replicates.
#' @param seed_base Integer seed offset; replicate \eqn{i} is seeded as
#'   \code{seed_base + i}.
#' @param parallel Logical; if \code{TRUE} and the \pkg{parallel}
#'   package is available, use \code{\link[parallel]{mclapply}}.
#'   Not used on Windows.
#' @param n_cores Number of worker processes when \code{parallel} is
#'   \code{TRUE}.
#' @param ... Additional arguments passed to the underlying simulator.
#' @return A list of length \code{n_reps}, each element the return
#'   value of \code{\link{simulate_trial}} or
#'   \code{\link{boin12_trial}}.
#' @examples
#' \dontrun{
#' scen <- ppod_scenarios$S4b
#' cell <- run_replicates(scen, design = "joint", n_reps = 20)
#' length(cell)
#' }
#' @export
run_replicates <- function(scen, design, n_reps = 400L, seed_base = 0L,
                            parallel = FALSE, n_cores = 2L, ...) {
  design <- match.arg(design, c("joint", "marginal", "plugin", "boin12"))
  one <- function(i) {
    set.seed(seed_base + i)
    if (design == "boin12") {
      boin12_trial(scen, ...)
    } else {
      simulate_trial(scen, criterion = design, ...)
    }
  }
  if (parallel &&
      requireNamespace("parallel", quietly = TRUE) &&
      .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(n_reps), one, mc.cores = n_cores)
  } else {
    lapply(seq_len(n_reps), one)
  }
}

#' Operating characteristics from a replicate cell
#'
#' Summarizes the output of \code{\link{run_replicates}} into the
#' operating characteristics used in the manuscript: probability of
#' correct OBD selection, mean utility loss against the true
#' \eqn{\Psi_j(\rho_j)} surface, probability of selecting a dose within
#' a fixed fraction of the utility range, probability of no selection,
#' allocation by exposure category (overdose, optimal, suboptimal),
#' and mean realized toxicity burden.
#'
#' @param cell Replicate list from \code{\link{run_replicates}}.
#' @param scen Scenario list (must have a \code{true_OBD} field; if
#'   absent, \code{true_OBD} is inferred from
#'   \code{\link{true_expected_utility}}).
#' @param eps Vector of utility-range fractions for near-optimal
#'   selection probability (default \code{c(0.05, 0.10)}).
#' @return A list with components:
#'   \describe{
#'     \item{\code{p_correct}, \code{se_correct}}{Estimated probability of
#'       correct OBD selection and its Monte Carlo standard error.}
#'     \item{\code{p_near}}{Named numeric vector of near-optimal
#'       selection probabilities at each \code{eps}.}
#'     \item{\code{p_no_select}}{Probability of no selection.}
#'     \item{\code{mean_loss}}{Mean utility loss against the true
#'       \eqn{\Psi_j(\rho_j)} surface. No-selection replicates
#'       contribute \eqn{\max_j \Psi_j(\rho_j)}.}
#'     \item{\code{selection_table}}{Length-\eqn{J + 1} vector of
#'       selection counts, with the last entry counting no-selection.}
#'     \item{\code{alloc}}{Named numeric vector \code{c(over, optimal,
#'       suboptimal)} of allocation proportions.}
#'     \item{\code{mean_tox_burden}}{Mean toxicity burden experienced
#'       by trial patients (working-scale toxicity outcomes for
#'       continuous designs, binary toxicity event rate for BOIN12).}
#'   }
#' @export
operating_characteristics <- function(cell, scen, eps = c(0.05, 0.10)) {
  J <- length(scen$mu_T)
  true_eu <- true_expected_utility(scen)
  true_obd <- if (!is.null(scen$true_OBD)) scen$true_OBD else which.max(true_eu)
  best_u <- true_eu[true_obd]
  util_range <- max(true_eu) - min(true_eu)

  # Tox burden: continuous designs have Y_T draws; BOIN12 has nT/n_dose.
  is_binary <- all(c("nT", "n_dose") %in% names(cell[[1L]]))

  sels <- vapply(cell, function(r) {
    if (is.null(r$selected) || is.na(r$selected)) NA_integer_ else r$selected
  }, integer(1))
  no_sel <- is.na(sels)

  losses <- numeric(length(cell))
  for (i in seq_along(cell)) {
    losses[i] <- if (no_sel[i]) best_u else best_u - true_eu[sels[i]]
  }

  p_correct <- mean(sels == true_obd, na.rm = TRUE) * mean(!no_sel)
  se_correct <- sqrt(p_correct * (1 - p_correct) / length(cell))

  p_near <- vapply(eps, function(e) {
    mean(!no_sel & losses <= e * util_range)
  }, numeric(1))
  names(p_near) <- paste0("eps_", eps)

  sel_table <- tabulate(sels, nbins = J)
  sel_table <- c(sel_table, no_select = sum(no_sel))

  # Allocation by exposure category.
  is_over <- scen$mu_T > 0                # working-scale overdose
  alloc_over <- alloc_opt <- alloc_sub <- numeric(length(cell))
  tox_burden <- numeric(length(cell))
  for (i in seq_along(cell)) {
    r <- cell[[i]]
    if (is_binary) {
      n <- sum(r$n_dose)
      alloc_over[i] <- sum(r$n_dose[is_over]) / max(n, 1L)
      alloc_opt[i]  <- r$n_dose[true_obd] / max(n, 1L)
      alloc_sub[i]  <- 1 - alloc_over[i] - alloc_opt[i]
      tox_burden[i] <- sum(r$nT) / max(n, 1L)
    } else {
      hist <- r$dose_history
      n <- length(hist)
      alloc_over[i] <- sum(is_over[hist]) / max(n, 1L)
      alloc_opt[i]  <- sum(hist == true_obd) / max(n, 1L)
      alloc_sub[i]  <- 1 - alloc_over[i] - alloc_opt[i]
      tox_burden[i] <- if (length(r$Y_T)) mean(r$Y_T) else NA_real_
    }
  }

  list(
    p_correct        = p_correct,
    se_correct       = se_correct,
    p_near           = p_near,
    p_no_select      = mean(no_sel),
    mean_loss        = mean(losses),
    selection_table  = sel_table,
    alloc            = c(over = mean(alloc_over),
                         optimal = mean(alloc_opt),
                         suboptimal = mean(alloc_sub)),
    mean_tox_burden  = mean(tox_burden, na.rm = TRUE)
  )
}
