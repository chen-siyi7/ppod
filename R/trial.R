#' Simulate one adaptive trial under a scenario
#'
#' Runs a single replicate of the adaptive design described in the
#' manuscript: an \eqn{n_0}-patient run-in phase using a posterior-mean
#' toxicity-exceedance criterion, followed by predictive-admissibility
#' gating and posterior predictive expected utility maximization on the
#' admissible set, with no-skipping-upward enforced throughout. The
#' trial stops early for excessive toxicity at \eqn{d_1} or for global
#' futility (empty admissible set).
#'
#' @param scen A scenario list, typically one element of
#'   \code{\link{ppod_scenarios}}, or built by \code{\link{make_scenario}}.
#'   Required fields are \code{mu_T}, \code{mu_E}, \code{rho},
#'   \code{sigma_T}, \code{sigma_E}, \code{util}, \code{pi_T},
#'   \code{pi_E}.
#' @param criterion Decision criterion: one of \code{"joint"},
#'   \code{"marginal"}, or \code{"plugin"}. Selects which expected-utility
#'   evaluator is used inside the dose-assignment rule and at final
#'   selection.
#' @param N_trial Total trial sample size (default 54).
#' @param cohort Cohort size (default 3).
#' @param n_runin Number of run-in patients before utility-based
#'   decisions begin (default 9).
#' @param S_interim,M_interim Posterior and predictive sample sizes for
#'   each interim decision (default 150 and 40).
#' @param S_final,M_final Posterior and predictive sample sizes at
#'   final OBD selection (default 400 and 100).
#' @param xT_max,xE_min Bounded-scale cutoffs used in the predictive
#'   admissibility gate (default 0.5 and \code{plogis(-1)}).
#' @param tau_T Working-scale toxicity cutoff for the run-in and safety
#'   stopping rules (default 0).
#' @param pi_T_stop Posterior probability cutoff for safety stopping at
#'   \eqn{d_1} (default 0.85).
#' @param verbose Logical, whether to print progress messages.
#' @return A list with elements \code{selected} (integer dose index, or
#'   \code{NA_integer_} if the trial stopped without a recommendation),
#'   \code{dose_history}, \code{Y_T}, \code{Y_E}, and logical flags
#'   \code{stopped_tox} and \code{stopped_fut}.
#' @examples
#' \donttest{
#' set.seed(1)
#' scen <- ppod_scenarios$S4b
#' res <- simulate_trial(scen, criterion = "joint")
#' res$selected
#' }
#' @export
simulate_trial <- function(scen, criterion = c("joint", "marginal", "plugin"),
                            N_trial = 54L, cohort = 3L, n_runin = 9L,
                            S_interim = 150L, M_interim = 40L,
                            S_final = 400L, M_final = 100L,
                            xT_max = 0.5, xE_min = stats::plogis(-1),
                            tau_T = 0, pi_T_stop = 0.85,
                            verbose = FALSE) {
  criterion <- match.arg(criterion)
  J <- length(scen$mu_T)
  rho_vec <- if (length(scen$rho) == 1L) rep(scen$rho, J) else scen$rho

  Y_T <- numeric(0)
  Y_E <- numeric(0)
  doses <- integer(0)
  dose_history <- integer(0)
  current <- 1L
  n_total <- 0L
  stopped_tox <- FALSE
  stopped_fut <- FALSE
  theta_hat <- NULL

  while (n_total < N_trial) {
    # Generate the next cohort at the current dose.
    for (k in seq_len(cohort)) {
      if (n_total >= N_trial) break
      z1 <- stats::rnorm(1)
      z2 <- stats::rnorm(1)
      yT <- scen$mu_T[current] + scen$sigma_T * z1
      yE <- scen$mu_E[current] +
            scen$sigma_E * (rho_vec[current] * z1 +
                            sqrt(max(1 - rho_vec[current]^2, 0)) * z2)
      Y_T <- c(Y_T, yT)
      Y_E <- c(Y_E, yE)
      doses <- c(doses, current)
      dose_history <- c(dose_history, current)
      n_total <- n_total + 1L
    }

    # Fit posterior.
    fit <- fit_laplace(Y_T, Y_E, doses, J = J, theta_init = theta_hat)
    if (!fit$ok) next
    theta_hat <- fit$theta
    samples <- sample_posterior(theta_hat, fit$cov, S_interim)

    # Safety stopping at d_1.
    post_excess <- .posterior_mean_tox_exceed(samples, tau_T, J)
    if (post_excess[1L] > pi_T_stop) {
      stopped_tox <- TRUE
      break
    }

    # Admissibility (only after the run-in; activity gate suspended at
    # very small per-dose sample sizes where the predictive activity
    # probability would be dominated by prior variance).
    if (n_total >= n_runin) {
      ad <- admissibility(samples, M_interim, xT_max, xE_min, J)
      n_per <- tabulate(doses, nbins = J)
      activity_ok <- (ad$p_under < scen$pi_E) | (n_per < 3L)
      admissible <- (ad$p_over < scen$pi_T) & activity_ok
      if (!any(admissible)) {
        stopped_fut <- TRUE
        break
      }
    } else {
      admissible <- rep(TRUE, J)
    }

    if (n_total >= N_trial) break

    if (n_total < n_runin) {
      # Run-in escalation: conservative posterior-mean rule.
      pe <- post_excess[current]
      if (pe < 0.10 && current < J) {
        current <- current + 1L
      } else if (pe > 0.40 && current > 1L) {
        current <- current - 1L
      }
    } else {
      # Utility-based assignment, with no-skipping-upward.
      highest <- max(doses)
      allowed <- admissible & (seq_len(J) <= highest + 1L)
      if (!any(allowed)) {
        stopped_fut <- TRUE
        break
      }
      eu <- switch(criterion,
        joint    = eu_joint(samples, M_interim, scen$util, J),
        marginal = eu_marginal(samples, M_interim, scen$util, J),
        plugin   = eu_plugin(theta_hat, scen$util, J)
      )
      eu_allowed <- ifelse(allowed, eu, -Inf)
      current <- as.integer(which.max(eu_allowed))
    }
    if (verbose) message("n = ", n_total, ", next dose = ", current)
  }

  # Final selection.
  if (stopped_tox || stopped_fut) {
    return(list(
      selected = NA_integer_,
      dose_history = dose_history, Y_T = Y_T, Y_E = Y_E,
      stopped_tox = stopped_tox, stopped_fut = stopped_fut
    ))
  }

  fit <- fit_laplace(Y_T, Y_E, doses, J = J, theta_init = theta_hat)
  if (!fit$ok) {
    return(list(
      selected = NA_integer_,
      dose_history = dose_history, Y_T = Y_T, Y_E = Y_E,
      stopped_tox = FALSE, stopped_fut = TRUE
    ))
  }
  samples <- sample_posterior(fit$theta, fit$cov, S_final)
  ad <- admissibility(samples, M_final, xT_max, xE_min, J)
  n_per <- tabulate(doses, nbins = J)
  activity_ok <- (ad$p_under < scen$pi_E) | (n_per < 3L)
  admissible <- (ad$p_over < scen$pi_T) & activity_ok
  if (!any(admissible)) {
    return(list(
      selected = NA_integer_,
      dose_history = dose_history, Y_T = Y_T, Y_E = Y_E,
      stopped_tox = FALSE, stopped_fut = TRUE
    ))
  }
  eu <- switch(criterion,
    joint    = eu_joint(samples, M_final, scen$util, J),
    marginal = eu_marginal(samples, M_final, scen$util, J),
    plugin   = eu_plugin(fit$theta, scen$util, J)
  )
  eu_a <- ifelse(admissible, eu, -Inf)
  list(
    selected = as.integer(which.max(eu_a)),
    dose_history = dose_history, Y_T = Y_T, Y_E = Y_E,
    stopped_tox = FALSE, stopped_fut = FALSE
  )
}
