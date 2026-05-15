#' Posterior predictive admissibility probabilities
#'
#' For each dose, returns the posterior predictive probability that a
#' future patient's bounded-scale toxicity exceeds a clinical cutoff
#' \eqn{x_T^{\max}}, and the posterior predictive probability that a
#' future patient's bounded-scale efficacy falls below a clinical
#' cutoff \eqn{x_E^{\min}}. These probabilities are computed from the
#' same Monte Carlo draws used by \code{\link{eu_joint}}, so the
#' admissibility gate and the utility evaluation share the same
#' posterior predictive sampling backbone.
#'
#' A dose is admissible if its over-toxicity probability is below
#' \eqn{\pi_T} and its under-efficacy probability is below
#' \eqn{\pi_E}; the exact gating rules are applied by
#' \code{\link{simulate_trial}}.
#'
#' @inheritParams eu_joint
#' @param xT_max Upper toxicity cutoff on the bounded scale.
#' @param xE_min Lower efficacy cutoff on the bounded scale.
#' @return A list with two length-\eqn{J} numeric vectors,
#'   \code{p_over} and \code{p_under}, giving the over-toxicity and
#'   under-efficacy posterior predictive probabilities at each dose.
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 30
#' doses <- rep(1:3, each = 10)
#' Y_T <- rnorm(n, mean = c(-1.5, -0.8, -0.2)[doses], sd = 0.5)
#' Y_E <- rnorm(n, mean = c(-0.5, 0.5, 0.7)[doses], sd = 0.5)
#' fit <- fit_laplace(Y_T, Y_E, doses, J = 3)
#' samples <- sample_posterior(fit$theta, fit$cov, S = 200)
#' admissibility(samples, M = 50, xT_max = 0.5,
#'               xE_min = plogis(-1), J = 3)
#' }
#' @export
admissibility <- function(samples, M, xT_max, xE_min, J) {
  up <- unpack_samples(samples, J)
  S <- nrow(samples)
  p_over <- numeric(J)
  p_under <- numeric(J)
  for (j in seq_len(J)) {
    z1 <- matrix(stats::rnorm(S * M), S, M)
    z2 <- matrix(stats::rnorm(S * M), S, M)
    s_omr <- sqrt(pmax(1 - up$rho^2, 0))
    yT <- up$mu_T[, j] + up$sT * z1
    yE <- up$mu_E[, j] + up$sE * (up$rho * z1 + s_omr * z2)
    xT <- stats::plogis(yT)
    xE <- stats::plogis(yE)
    p_over[j]  <- mean(xT > xT_max)
    p_under[j] <- mean(xE < xE_min)
  }
  list(p_over = p_over, p_under = p_under)
}

# Posterior probability that mu_T(d_j) exceeds a working-scale toxicity
# cutoff, marginalized over posterior parameter uncertainty. Used by
# the run-in escalation rule and the safety stopping rule, where the
# predictive distribution would be uninformative at very small per-dose
# sample sizes.
.posterior_mean_tox_exceed <- function(samples, tau_T, J) {
  up <- unpack_samples(samples, J)
  colMeans(up$mu_T > tau_T)
}
