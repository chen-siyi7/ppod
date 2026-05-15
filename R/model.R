# Model parameterization and log-posterior for the joint dose-response
# model of Section 2 of the manuscript.

#' Unpack a parameter vector for the joint dose-response model
#'
#' The joint Gaussian dose-response model has \eqn{2J + 3} parameters
#' organized as
#' \deqn{\theta = (\alpha_T, \beta_{T,2}, \ldots, \beta_{T,J}, \alpha_E,
#' \tilde\gamma_1, \ldots, \tilde\gamma_{J-1},
#' \log\sigma_T, \log\sigma_E, \mathrm{atanh}\,\rho).}
#' This function reverses the packing and returns the components on
#' interpretable scales, including \eqn{\gamma_j} on the sum-to-zero
#' constraint (the last component is \eqn{-\sum_{j=1}^{J-1} \tilde\gamma_j}).
#'
#' @param theta Numeric vector of length \eqn{2J + 3}.
#' @param J Integer, number of dose levels.
#' @return A list with components \code{alpha_T}, \code{beta_T} (length
#'   \eqn{J-1}), \code{alpha_E}, \code{gamma} (length \eqn{J},
#'   sum-to-zero), \code{log_sT}, \code{log_sE}, \code{atanh_rho},
#'   \code{sT}, \code{sE}, and \code{rho}.
#' @examples
#' theta <- rep(0, 13)              # 2*5 + 3 = 13 parameters for J = 5
#' unpack_theta(theta, J = 5)
#' @export
unpack_theta <- function(theta, J) {
  stopifnot(length(theta) == 2L * J + 3L, J >= 3L)
  alpha_T   <- theta[1L]
  beta_T    <- theta[2L:J]                          # length J - 1
  alpha_E   <- theta[J + 1L]
  gamma_fr  <- theta[(J + 2L):(2L * J)]             # length J - 1
  gamma     <- c(gamma_fr, -sum(gamma_fr))          # sum-to-zero
  log_sT    <- theta[2L * J + 1L]
  log_sE    <- theta[2L * J + 2L]
  atanh_rho <- theta[2L * J + 3L]
  list(
    alpha_T   = alpha_T,
    beta_T    = beta_T,
    alpha_E   = alpha_E,
    gamma     = gamma,
    log_sT    = log_sT,
    log_sE    = log_sE,
    atanh_rho = atanh_rho,
    sT        = exp(log_sT),
    sE        = exp(log_sE),
    rho       = tanh(atanh_rho)
  )
}

# Build mu_T(d_1), ..., mu_T(d_J) from alpha_T and beta_T using the
# softplus-increment parameterization that enforces monotonicity.
.mu_T_vec <- function(alpha_T, beta_T) {
  alpha_T + c(0, cumsum(softplus(beta_T)))
}

#' Negative log posterior for the joint dose-response model
#'
#' Computes the negative joint log-posterior at parameter vector
#' \code{theta} given paired toxicity and efficacy observations at known
#' doses. The likelihood is bivariate Gaussian with dose-level means
#' built from a softplus-increment monotone toxicity model and a
#' second-order random-walk efficacy model under a sum-to-zero
#' identifiability constraint. Priors are described in
#' \code{vignette("ppod-overview")}.
#'
#' This is intended primarily as the objective passed to
#' \code{\link[stats]{optim}} inside \code{\link{fit_laplace}}.
#'
#' @param theta Numeric parameter vector of length \eqn{2J + 3}.
#' @param Y_T,Y_E Numeric vectors of observed toxicity and efficacy on
#'   the working scale, same length \eqn{n}.
#' @param doses Integer vector of dose indices (1-based), length \eqn{n}.
#' @param J Integer, number of dose levels.
#' @return A single non-negative number, the negative log posterior up
#'   to a constant. Returns a large finite value if the proposed
#'   parameters violate \eqn{|\rho| < 1}.
#' @examples
#' set.seed(1)
#' Y_T <- rnorm(9, mean = c(-1.5, -0.5, 0.5)[rep(1:3, each = 3)])
#' Y_E <- rnorm(9, mean = c(-0.5, 0.5, 1.0)[rep(1:3, each = 3)])
#' doses <- rep(1:3, each = 3)
#' theta <- rep(0, 9)               # 2*3 + 3 = 9
#' neg_log_post(theta, Y_T, Y_E, doses, J = 3)
#' @export
neg_log_post <- function(theta, Y_T, Y_E, doses, J) {
  p <- unpack_theta(theta, J)
  mu_T <- .mu_T_vec(p$alpha_T, p$beta_T)
  mu_E <- p$alpha_E + p$gamma
  one_m <- 1 - p$rho^2
  if (one_m <= 1e-12) {
    return(1e20)
  }
  eT <- Y_T - mu_T[doses]
  eE <- Y_E - mu_E[doses]
  n  <- length(Y_T)
  log_det <- log(p$sT^2 * p$sE^2 * one_m)
  quad <- (eT^2 / p$sT^2
           - 2 * p$rho * eT * eE / (p$sT * p$sE)
           + eE^2 / p$sE^2) / one_m
  log_lik <- -0.5 * n * (2 * log(2 * pi) + log_det) - 0.5 * sum(quad)

  # Priors
  lp <- -0.5 * p$alpha_T^2 / 4                                  # N(0, 4)
  lp <- lp - 0.5 * sum((p$beta_T - (-1))^2 / 2)                 # N(-1, 2)
  lp <- lp - 0.5 * p$alpha_E^2 / 4                              # N(0, 4)
  if (J >= 3L) {
    d2 <- p$gamma[3L:J] - 2 * p$gamma[2L:(J - 1L)] + p$gamma[1L:(J - 2L)]
    lp <- lp - 0.5 * sum(d2^2)                                  # RW2, tau = 1
  }
  lp <- lp - 0.5 * p$gamma[1L]^2 / 4 - 0.5 * p$gamma[2L]^2 / 4  # boundary
  lp <- lp + p$log_sT - 0.5 * exp(2 * p$log_sT)                 # Half-N(0,1)
  lp <- lp + p$log_sE - 0.5 * exp(2 * p$log_sE)
  lp <- lp + log(one_m)                                         # atanh Jac.

  -(log_lik + lp)
}

# Numerical Hessian using central differences for diagonal entries and
# the four-point cross formula for off-diagonal entries. Called only
# from fit_laplace at the posterior mode; fast at p <= 13.
.numerical_hessian <- function(f, x, args, eps = 1e-4) {
  n <- length(x)
  f0 <- do.call(f, c(list(x), args))
  fp <- fm <- numeric(n)
  for (i in seq_len(n)) {
    xi <- x; xi[i] <- xi[i] + eps; fp[i] <- do.call(f, c(list(xi), args))
    xi <- x; xi[i] <- xi[i] - eps; fm[i] <- do.call(f, c(list(xi), args))
  }
  H <- matrix(0, n, n)
  for (i in seq_len(n)) H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / eps^2
  if (n >= 2L) {
    for (i in seq_len(n - 1L)) {
      for (j in (i + 1L):n) {
        xij <- x; xij[i] <- xij[i] + eps; xij[j] <- xij[j] + eps
        fpp <- do.call(f, c(list(xij), args))
        H[i, j] <- (fpp - fp[i] - fp[j] + f0) / eps^2
        H[j, i] <- H[i, j]
      }
    }
  }
  H
}
