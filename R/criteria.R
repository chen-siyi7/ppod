#' Unpack posterior samples into mean and dispersion components
#'
#' Given a matrix of posterior samples of \eqn{\theta}, returns the
#' implied dose-level toxicity and efficacy means and the within-dose
#' standard deviations and correlation, on interpretable scales.
#' Used internally by the predictive utility and admissibility
#' evaluators but exposed for users implementing custom criteria.
#'
#' @param samples \eqn{S}-by-\eqn{2J + 3} matrix of posterior samples,
#'   typically returned by \code{\link{sample_posterior}}.
#' @param J Integer, number of dose levels.
#' @return A list with components \code{mu_T} and \code{mu_E}, each an
#'   \eqn{S}-by-\eqn{J} matrix of dose-level mean draws, and length-\eqn{S}
#'   vectors \code{sT}, \code{sE}, and \code{rho}.
#' @export
unpack_samples <- function(samples, J) {
  S <- nrow(samples)
  aT     <- samples[, 1L]
  bT     <- samples[, 2L:J, drop = FALSE]
  aE     <- samples[, J + 1L]
  gf     <- samples[, (J + 2L):(2L * J), drop = FALSE]
  gamma  <- cbind(gf, -rowSums(gf))                  # S x J, sum-to-zero
  log_sT <- samples[, 2L * J + 1L]
  log_sE <- samples[, 2L * J + 2L]
  atRho  <- samples[, 2L * J + 3L]
  incs   <- log1p(exp(bT))                           # S x (J-1)
  cum    <- t(apply(incs, 1L, cumsum))               # S x (J-1)
  if (ncol(cum) == 1L) cum <- matrix(cum, S, 1L)
  mu_T   <- cbind(aT, aT + cum)                      # S x J
  mu_E   <- aE + gamma                               # S x J
  list(mu_T = mu_T, mu_E = mu_E,
       sT = exp(log_sT), sE = exp(log_sE), rho = tanh(atRho))
}

#' Posterior predictive expected utility under the joint distribution
#'
#' Monte Carlo estimate of \eqn{\widehat U_n(d_j) =
#' E[U(g_T(Y_T), g_E(Y_E)) \mid d_j, \mathcal D_n]} at each dose
#' \eqn{d_j}, computed under the joint Gaussian working model. For each
#' posterior draw \eqn{\theta^{(s)}}, \code{M} predictive outcome pairs
#' are drawn from \eqn{\mathcal N_2(\mu^{(s)}(d_j), \Sigma_j^{(s)})}
#' with the within-dose correlation \eqn{\rho^{(s)}} applied; the
#' utility is averaged over the \eqn{SM} draws. This is the joint
#' posterior predictive criterion of the manuscript.
#'
#' @param samples Posterior sample matrix from \code{\link{sample_posterior}}.
#' @param M Number of predictive outcome draws per posterior sample.
#' @param util Utility function returning \eqn{U(x_T, x_E)} on the bounded
#'   scale, typically built by \code{\link{make_utility}}.
#' @param J Integer, number of dose levels.
#' @return Numeric vector of length \eqn{J} giving estimated expected
#'   utility at each dose.
#' @seealso \code{\link{eu_marginal}}, \code{\link{eu_plugin}}.
#' @export
eu_joint <- function(samples, M, util, J) {
  up <- unpack_samples(samples, J)
  S <- nrow(samples)
  eu <- numeric(J)
  for (j in seq_len(J)) {
    z1 <- matrix(stats::rnorm(S * M), S, M)
    z2 <- matrix(stats::rnorm(S * M), S, M)
    s_omr <- sqrt(pmax(1 - up$rho^2, 0))
    yT <- up$mu_T[, j] + up$sT * z1
    yE <- up$mu_E[, j] + up$sE * (up$rho * z1 + s_omr * z2)
    xT <- stats::plogis(yT)
    xE <- stats::plogis(yE)
    eu[j] <- mean(util(xT, xE))
  }
  eu
}

#' Posterior predictive expected utility under the marginal distribution
#'
#' Monte Carlo estimate of expected utility at each dose, evaluated as
#' if \eqn{Y_T} and \eqn{Y_E} were independent within dose. The
#' underlying dose-response model is the same joint Gaussian
#' specification as for \code{\link{eu_joint}}; the only difference is
#' that the Monte Carlo draws set \eqn{\rho} to zero at the utility
#' integration step. The contrast between \code{eu_joint} and
#' \code{eu_marginal} isolates the role of dependence-sensitive
#' integration in the posterior predictive criterion.
#'
#' @inheritParams eu_joint
#' @return Numeric vector of length \eqn{J}.
#' @seealso \code{\link{eu_joint}}, \code{\link{eu_plugin}}.
#' @export
eu_marginal <- function(samples, M, util, J) {
  up <- unpack_samples(samples, J)
  S <- nrow(samples)
  eu <- numeric(J)
  for (j in seq_len(J)) {
    z1 <- matrix(stats::rnorm(S * M), S, M)
    z2 <- matrix(stats::rnorm(S * M), S, M)
    yT <- up$mu_T[, j] + up$sT * z1
    yE <- up$mu_E[, j] + up$sE * z2
    xT <- stats::plogis(yT)
    xE <- stats::plogis(yE)
    eu[j] <- mean(util(xT, xE))
  }
  eu
}

#' Plug-in expected utility at the posterior mean
#'
#' Evaluates \eqn{U(g_T(\hat\mu_T(d_j)), g_E(\hat\mu_E(d_j)))} at the
#' posterior mode \eqn{\hat\theta}. This plug-in criterion ignores both
#' residual outcome variability and posterior parameter uncertainty,
#' which makes it a useful comparator for isolating the role of
#' posterior predictive integration in OBD selection.
#'
#' @param theta_hat Posterior mode from \code{\link{fit_laplace}}, a
#'   vector of length \eqn{2J + 3}.
#' @param util Utility function returning \eqn{U(x_T, x_E)} on the
#'   bounded scale.
#' @param J Integer, number of dose levels.
#' @return Numeric vector of length \eqn{J}.
#' @seealso \code{\link{eu_joint}}, \code{\link{eu_marginal}}.
#' @export
eu_plugin <- function(theta_hat, util, J) {
  p <- unpack_theta(theta_hat, J)
  mu_T <- .mu_T_vec(p$alpha_T, p$beta_T)
  mu_E <- p$alpha_E + p$gamma
  xT <- stats::plogis(mu_T)
  xE <- stats::plogis(mu_E)
  util(xT, xE)
}

#' Draw posterior predictive outcomes at a fixed dose
#'
#' Convenience function for downstream users who want raw posterior
#' predictive outcome draws (not just an expected utility). Returns
#' \eqn{(Y_T, Y_E)} pairs on the working scale and \eqn{(x_T, x_E)}
#' pairs on the bounded scale at a specified dose.
#'
#' @inheritParams eu_joint
#' @param dose Integer dose index in \code{1:J}.
#' @return A list with matrices \code{Y_T}, \code{Y_E}, \code{x_T},
#'   \code{x_E}, each \eqn{S}-by-\eqn{M}.
#' @export
predict_outcomes <- function(samples, M, dose, J) {
  stopifnot(dose >= 1L, dose <= J)
  up <- unpack_samples(samples, J)
  S <- nrow(samples)
  z1 <- matrix(stats::rnorm(S * M), S, M)
  z2 <- matrix(stats::rnorm(S * M), S, M)
  s_omr <- sqrt(pmax(1 - up$rho^2, 0))
  yT <- up$mu_T[, dose] + up$sT * z1
  yE <- up$mu_E[, dose] + up$sE * (up$rho * z1 + s_omr * z2)
  list(
    Y_T = yT, Y_E = yE,
    x_T = stats::plogis(yT), x_E = stats::plogis(yE)
  )
}
