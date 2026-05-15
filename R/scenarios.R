#' Construct a scenario list
#'
#' Helper to build a scenario list in the format expected by
#' \code{\link{simulate_trial}}, \code{\link{boin12_trial}}, and
#' \code{\link{true_expected_utility}}. Users typically use one of the
#' built-in scenarios in \code{\link{ppod_scenarios}}; this constructor
#' is for users who want to define their own.
#'
#' @param mu_T Numeric vector of working-scale toxicity means, length
#'   \eqn{J}.
#' @param mu_E Numeric vector of working-scale efficacy means, length
#'   \eqn{J}.
#' @param rho Either a single correlation (applied to all doses) or a
#'   length-\eqn{J} vector of dose-specific correlations.
#' @param sigma_T,sigma_E Working-scale residual standard deviations.
#' @param util Utility function on the bounded scale, typically built
#'   by \code{\link{make_utility}}.
#' @param pi_T,pi_E Tail-probability cutoffs for the predictive
#'   admissibility gate (default 0.20 and 0.30).
#' @param true_OBD Optional integer giving the dose label for the true
#'   OBD; if \code{NULL}, computed by \code{\link{true_expected_utility}}.
#' @return A scenario list.
#' @examples
#' s <- make_scenario(
#'   mu_T = c(-2.2, -1.5, -0.8, -0.2, 0.4),
#'   mu_E = c(-1.4, -0.2, 0.6, 1.4, 2.2),
#'   rho  = 0
#' )
#' s$true_OBD
#' @export
make_scenario <- function(mu_T, mu_E, rho = 0,
                           sigma_T = 0.5, sigma_E = 0.5,
                           util = make_utility(1, 1, 2),
                           pi_T = 0.20, pi_E = 0.30,
                           true_OBD = NULL) {
  stopifnot(length(mu_T) == length(mu_E))
  stopifnot(length(rho) == 1L || length(rho) == length(mu_T))
  scen <- list(
    mu_T = mu_T, mu_E = mu_E,
    rho = rho, sigma_T = sigma_T, sigma_E = sigma_E,
    util = util, pi_T = pi_T, pi_E = pi_E
  )
  if (is.null(true_OBD)) {
    eu <- true_expected_utility(scen)
    scen$true_OBD <- as.integer(which.max(eu))
  } else {
    scen$true_OBD <- as.integer(true_OBD)
  }
  scen
}

#' Simulation scenarios used in the manuscript
#'
#' A list of 12 scenarios reproducing the simulation study in the
#' accompanying manuscript. Each element is a scenario list of the form
#' returned by \code{\link{make_scenario}}.
#'
#' Defined programmatically in \code{R/scenarios-data.R} and built when
#' the package is loaded.
#'
#' @format A named list of scenarios:
#' \describe{
#'   \item{S1}{Monotone efficacy.}
#'   \item{S2}{Plateau efficacy.}
#'   \item{S3}{Unimodal efficacy with mild downturn.}
#'   \item{S4a, S4b, S4c}{Unimodal efficacy at \eqn{\rho \in \{-0.5, 0, 0.5\}};
#'     dependence-sensitivity contrast under constant \eqn{\rho}.}
#'   \item{S4d}{Dose-varying dependence stress test:
#'     \eqn{\rho = (-0.85, -0.85, 0.95, 0.92, 0.9)}, narrow marginal gap
#'     between \eqn{d_2} and \eqn{d_3}, larger residual scale. The
#'     deterministic joint OBD is \eqn{d_2} versus marginal OBD \eqn{d_3}.}
#'   \item{S5}{Unimodal immuno-oncology profile.}
#'   \item{S6}{Higher-toxicity targeted-agent profile.}
#'   \item{S7a, S7b, S7c}{Close-tied marginal preference at
#'     \eqn{\rho \in \{-0.5, 0, 0.5\}}; deterministic joint OBD is
#'     \eqn{d_3} under negative dependence and \eqn{d_2} otherwise.}
#' }
#' @source Manuscript simulation specifications.
#' @examples
#' names(ppod_scenarios)
#' ppod_scenarios$S4d$rho
#' ppod_scenarios$S4d$true_OBD
#' @name ppod_scenarios
NULL
