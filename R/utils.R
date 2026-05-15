#' Softplus function
#'
#' Smooth, non-negative function \code{log(1 + exp(x))}. Used to
#' parameterize non-negative dose increments in the monotone toxicity
#' model: each dose-level increment in the toxicity mean is a softplus
#' transformation of an unconstrained real parameter.
#'
#' Implemented in a numerically stable way via \code{\link{log1p}}.
#'
#' @param x Numeric vector.
#' @return A numeric vector of the same length as \code{x} with
#'   non-negative values.
#' @examples
#' softplus(c(-5, 0, 5))
#' @export
softplus <- function(x) {
  log1p(exp(x))
}

#' Logistic (sigmoid) transformation
#'
#' Convenience wrapper for \code{\link[stats]{plogis}}. Used as the
#' bounded transformation \eqn{g_T = g_E = \mathrm{logit}^{-1}} that maps
#' the unbounded working-scale outcomes to the open unit interval on
#' which the utility function is defined.
#'
#' @param y Numeric vector on the working scale.
#' @return Numeric vector on \code{(0, 1)}.
#' @examples
#' g_logistic(c(-2, 0, 2))
#' @export
g_logistic <- function(y) {
  stats::plogis(y)
}

# Number of model parameters given J doses.
# alpha_T, beta_{T,2..J}, alpha_E, gamma_{1..J-1}, log sigma_T, log sigma_E,
# atanh rho. That is 2J + 3 parameters.
.n_params <- function(J) {
  2L * J + 3L
}

# Default theta initializer. Used by fit_laplace when no warm start is given.
.theta_default <- function(J) {
  P <- .n_params(J)
  theta <- numeric(P)
  theta[1L]            <- -1.5            # alpha_T
  theta[2L:J]          <- -1.0            # beta_{T,k}
  theta[J + 1L]        <-  0              # alpha_E
  # gamma_{1..J-1} stay at 0
  theta[2L * J + 1L]   <- log(0.5)        # log sigma_T
  theta[2L * J + 2L]   <- log(0.5)        # log sigma_E
  theta[2L * J + 3L]   <-  0              # atanh rho
  theta
}
