#' Deterministic true expected utility under a scenario
#'
#' Computes the true expected utility \eqn{\Psi_j(\rho_j) = E_{\theta^*}\{
#' U(g_T(Y_T), g_E(Y_E)) \mid d_j\}} at each dose under the
#' data-generating parameters of a scenario, by Monte Carlo over the
#' joint bivariate normal data-generating distribution. The
#' dose-maximizing \eqn{\Psi_j(\rho_j)} on the predictive admissible
#' set defines the true OBD; the marginal counterpart with
#' \eqn{\rho_j = 0} defines the marginal OBD. These two true values
#' are reported in the manuscript's Table 1 and underlie the
#' deterministic truth tables for the S4d and S7 scenarios.
#'
#' @param scen Scenario list.
#' @param n_mc Number of Monte Carlo draws per dose (default
#'   \eqn{2 \times 10^5}).
#' @param marginal Logical; if \code{TRUE}, evaluate under independence
#'   (\eqn{\rho_j \equiv 0}) instead of the scenario's \eqn{\rho_j}.
#'   This is useful for computing the marginal-OBD truth value.
#' @param seed Optional integer seed. If \code{NULL}, a stable hash of
#'   the scenario means is used so that calls on the same scenario
#'   return the same Monte Carlo estimate.
#' @return Numeric vector of length \code{length(scen$mu_T)} giving
#'   estimated expected utility at each dose.
#' @examples
#' true_expected_utility(ppod_scenarios$S4b)
#' true_expected_utility(ppod_scenarios$S4d)
#' true_expected_utility(ppod_scenarios$S4d, marginal = TRUE)
#' @export
true_expected_utility <- function(scen, n_mc = 200000L,
                                    marginal = FALSE, seed = NULL) {
  J <- length(scen$mu_T)
  rho_vec <- if (length(scen$rho) == 1L) rep(scen$rho, J) else scen$rho
  if (marginal) rho_vec <- rep(0, J)
  if (is.null(seed)) {
    seed <- 99L + as.integer(sum(scen$mu_T * 1000))
  }
  set.seed(seed)
  eu <- numeric(J)
  for (j in seq_len(J)) {
    z1 <- stats::rnorm(n_mc)
    z2 <- stats::rnorm(n_mc)
    yT <- scen$mu_T[j] + scen$sigma_T * z1
    yE <- scen$mu_E[j] +
          scen$sigma_E * (rho_vec[j] * z1 +
                          sqrt(max(1 - rho_vec[j]^2, 0)) * z2)
    eu[j] <- mean(scen$util(stats::plogis(yT), stats::plogis(yE)))
  }
  eu
}
