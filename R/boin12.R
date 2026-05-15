#' BOIN12 binary-endpoint comparator
#'
#' Simulates one trial under a simplified BOIN12 design (Lin et al.,
#' 2020) using binary toxicity \eqn{I(Y_T > \tau_T)} and binary
#' efficacy \eqn{I(Y_E > \tau_E)} derived by thresholding the
#' scenario's continuous outcomes. The binary cutoffs match the
#' clinically meaningful thresholds used by the proposed design's
#' predictive admissibility gate (no scenario-specific retuning). This
#' is the binary-endpoint comparator reported in the simulation results
#' of the manuscript; the continuous versus BOIN12 contrast quantifies
#' the combined effect of endpoint type, decision rule, efficacy model,
#' and admissibility mechanics.
#'
#' @param scen Scenario list, as for \code{\link{simulate_trial}}.
#' @param target_T,target_E BOIN12 target binary probabilities.
#' @param w11,w10,w01,w00 BOIN12 utility weights for the four
#'   (toxicity, efficacy) outcome configurations.
#' @param N_trial Total trial sample size (default 54).
#' @param cohort Cohort size (default 3).
#' @param tau_T,tau_E Working-scale thresholds defining the induced
#'   binary endpoints (default \eqn{0} and \eqn{-1}).
#' @return A list with components \code{selected}, \code{dose_history},
#'   \code{n_dose}, \code{nT}, \code{nE}, \code{stopped_tox},
#'   \code{stopped_fut}, mirroring the return value of
#'   \code{\link{simulate_trial}}.
#' @references
#' Lin, R., Zhou, Y., Yan, F., Li, D., and Yuan, Y. (2020). BOIN12:
#' Bayesian optimal interval phase I/II trial design for utility-based
#' dose finding in immunotherapy and targeted therapies. \emph{JCO
#' Precision Oncology}, 4, 1393-1402.
#' @export
boin12_trial <- function(scen,
                          target_T = 0.30, target_E = 0.40,
                          w11 = 1.0, w10 = 0.4, w01 = 0.3, w00 = 0.0,
                          N_trial = 54L, cohort = 3L,
                          tau_T = 0, tau_E = -1) {
  J <- length(scen$mu_T)
  rho_vec <- if (length(scen$rho) == 1L) rep(scen$rho, J) else scen$rho
  # Induced binary probabilities under the scenario's continuous model.
  pT_true <- 1 - stats::pnorm(tau_T, mean = scen$mu_T, sd = scen$sigma_T)
  pE_true <- 1 - stats::pnorm(tau_E, mean = scen$mu_E, sd = scen$sigma_E)

  nT <- nE <- n_dose <- integer(J)
  dose_history <- integer(0)
  current <- 1L
  n_total <- 0L

  esc_cut    <- 0.6 * target_T
  de_esc_cut <- 1.4 * target_T

  while (n_total < N_trial) {
    for (k in seq_len(cohort)) {
      if (n_total >= N_trial) break
      z <- MASS::mvrnorm(
        1L, mu = c(0, 0),
        Sigma = matrix(c(1, rho_vec[current], rho_vec[current], 1), 2L, 2L)
      )
      yT <- as.integer(stats::pnorm(z[1L]) < pT_true[current])
      yE <- as.integer(stats::pnorm(z[2L]) < pE_true[current])
      nT[current]      <- nT[current] + yT
      nE[current]      <- nE[current] + yE
      n_dose[current]  <- n_dose[current] + 1L
      dose_history     <- c(dose_history, current)
      n_total          <- n_total + 1L
    }

    # Safety stopping at d_1.
    if (n_dose[1L] > 0L) {
      a <- 1 + nT[1L]
      b <- 1 + n_dose[1L] - nT[1L]
      if (stats::pbeta(target_T, a, b, lower.tail = FALSE) > 0.95) {
        return(list(
          selected = NA_integer_, dose_history = dose_history,
          n_dose = n_dose, nT = nT, nE = nE,
          stopped_tox = TRUE, stopped_fut = FALSE
        ))
      }
    }

    if (n_total >= N_trial) break

    # BOIN12 interval-based dose move.
    pT_hat <- nT[current] / n_dose[current]
    if (pT_hat <= esc_cut && current < J) {
      if (n_dose[current + 1L] > 0L) {
        a <- 1 + nT[current + 1L]
        b <- 1 + n_dose[current + 1L] - nT[current + 1L]
        if (stats::pbeta(target_T, a, b, lower.tail = FALSE) <= 0.85) {
          current <- current + 1L
        }
      } else {
        current <- current + 1L
      }
    } else if (pT_hat >= de_esc_cut && current > 1L) {
      current <- current - 1L
    }
  }

  # End-of-trial admissibility and utility selection.
  admiss <- logical(J)
  util <- rep(-Inf, J)
  for (j in seq_len(J)) {
    if (n_dose[j] == 0L) next
    aT <- 1 + nT[j]; bT <- 1 + n_dose[j] - nT[j]
    aE <- 1 + nE[j]; bE <- 1 + n_dose[j] - nE[j]
    ok_T <- stats::pbeta(target_T, aT, bT, lower.tail = FALSE) < 0.85
    ok_E <- stats::pbeta(target_E, aE, bE) < 0.85
    admiss[j] <- ok_T && ok_E
    if (admiss[j]) {
      pT_p <- aT / (aT + bT)
      pE_p <- aE / (aE + bE)
      util[j] <- (1 - pT_p) * pE_p * w11 +
                 (1 - pT_p) * (1 - pE_p) * w10 +
                 pT_p * pE_p * w01 +
                 pT_p * (1 - pE_p) * w00
    }
  }
  if (!any(admiss)) {
    return(list(
      selected = NA_integer_, dose_history = dose_history,
      n_dose = n_dose, nT = nT, nE = nE,
      stopped_tox = FALSE, stopped_fut = TRUE
    ))
  }
  list(
    selected = as.integer(which.max(util)),
    dose_history = dose_history,
    n_dose = n_dose, nT = nT, nE = nE,
    stopped_tox = FALSE, stopped_fut = FALSE
  )
}
