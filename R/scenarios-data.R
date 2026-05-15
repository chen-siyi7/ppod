# Definition of ppod_scenarios. Each scenario carries its own utility
# function (closure) so we cannot store as .rda without serializing
# closures; defining in package source keeps the object reproducible
# and inspectable.

# These S4 / S7 utility families correspond to the manuscript:
#   S1-S6: (alpha, beta, lam) = (1, 1, 2)
#   S7   : (alpha, beta, lam) = (1, 2, 3)
.util_S4 <- make_utility(1, 1, 2)
.util_S7 <- make_utility(1, 2, 3)

#' @rdname ppod_scenarios
#' @export
ppod_scenarios <- list(
  S1   = list(
    mu_T = c(-2.2, -1.5, -0.8, -0.2, 0.4),
    mu_E = c(-1.4, -0.2,  0.6,  1.4, 2.2),
    rho = 0, sigma_T = 0.5, sigma_E = 0.5,
    util = .util_S4, pi_T = 0.20, pi_E = 0.30, true_OBD = 3L
  ),
  S2   = list(
    mu_T = c(-2.2, -1.5, -0.8, -0.2, 0.4),
    mu_E = c(-1.4, -0.2,  0.8,  0.95, 0.95),
    rho = 0, sigma_T = 0.5, sigma_E = 0.5,
    util = .util_S4, pi_T = 0.20, pi_E = 0.30, true_OBD = 3L
  ),
  S3   = list(
    mu_T = c(-2.2, -1.5, -0.8, -0.2, 0.4),
    mu_E = c(-1.4, -0.2,  0.8,  0.7, 0.2),
    rho = 0, sigma_T = 0.5, sigma_E = 0.5,
    util = .util_S4, pi_T = 0.20, pi_E = 0.30, true_OBD = 3L
  ),
  S4a  = list(
    mu_T = c(-2.2, -1.5, -0.8, -0.2, 0.4),
    mu_E = c(-1.4, -0.2,  0.8,  0.7, 0.2),
    rho = -0.5, sigma_T = 0.5, sigma_E = 0.5,
    util = .util_S4, pi_T = 0.20, pi_E = 0.30, true_OBD = 3L
  ),
  S4b  = list(
    mu_T = c(-2.2, -1.5, -0.8, -0.2, 0.4),
    mu_E = c(-1.4, -0.2,  0.8,  0.7, 0.2),
    rho = 0, sigma_T = 0.5, sigma_E = 0.5,
    util = .util_S4, pi_T = 0.20, pi_E = 0.30, true_OBD = 3L
  ),
  S4c  = list(
    mu_T = c(-2.2, -1.5, -0.8, -0.2, 0.4),
    mu_E = c(-1.4, -0.2,  0.8,  0.7, 0.2),
    rho = 0.5, sigma_T = 0.5, sigma_E = 0.5,
    util = .util_S4, pi_T = 0.20, pi_E = 0.30, true_OBD = 3L
  ),
  S4d  = list(
    mu_T = c(-2.2, -1.5, -0.8, -0.2, 0.4),
    mu_E = c(-1.4, -0.1,  0.7,  0.42, -0.3),
    rho = c(-0.85, -0.85, 0.95, 0.92, 0.9),
    sigma_T = 0.8, sigma_E = 0.8,
    util = .util_S4, pi_T = 0.25, pi_E = 0.30, true_OBD = 2L
  ),
  S5   = list(
    mu_T = c(-2.9, -2.0, -1.1, -0.4, 0.3),
    mu_E = c(-1.4,  0.0,  1.1,  0.5, -0.2),
    rho = 0.3, sigma_T = 0.5, sigma_E = 0.5,
    util = .util_S4, pi_T = 0.20, pi_E = 0.30, true_OBD = 3L
  ),
  S6   = list(
    mu_T = c(-1.4, -0.8, -0.3,  0.2, 0.85),
    mu_E = c(-0.85, 0.2,  0.6,  0.5, 0.2),
    rho = 0.5, sigma_T = 0.5, sigma_E = 0.5,
    util = .util_S4, pi_T = 0.20, pi_E = 0.30, true_OBD = 2L
  ),
  S7a  = list(
    mu_T = c(-2.2, -1.0, -0.4,  0.2, 0.85),
    mu_E = c(-1.4,  0.0,  0.6,  0.3, -0.2),
    rho = -0.5, sigma_T = 0.7, sigma_E = 0.7,
    util = .util_S7, pi_T = 0.30, pi_E = 0.30, true_OBD = 3L
  ),
  S7b  = list(
    mu_T = c(-2.2, -1.0, -0.4,  0.2, 0.85),
    mu_E = c(-1.4,  0.0,  0.6,  0.3, -0.2),
    rho = 0, sigma_T = 0.7, sigma_E = 0.7,
    util = .util_S7, pi_T = 0.30, pi_E = 0.30, true_OBD = 2L
  ),
  S7c  = list(
    mu_T = c(-2.2, -1.0, -0.4,  0.2, 0.85),
    mu_E = c(-1.4,  0.0,  0.6,  0.3, -0.2),
    rho = 0.5, sigma_T = 0.7, sigma_E = 0.7,
    util = .util_S7, pi_T = 0.30, pi_E = 0.30, true_OBD = 2L
  )
)
