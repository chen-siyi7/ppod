#' Construct a substitution utility function
#'
#' Returns a utility function on the bounded unit square \eqn{(0,1)^2}
#' from the parametric family
#' \deqn{U(x_T, x_E) = x_E^\alpha \cdot \exp(-\lambda x_T^\beta),}
#' where \eqn{x_T} and \eqn{x_E} are toxicity and efficacy on the
#' bounded transformed scale. Under \eqn{\alpha, \beta \ge 1} and
#' \eqn{\lambda > 0}, the mixed partial \eqn{U_{TE}} is strictly
#' negative on the open unit square, which is the substitutability
#' condition under which Corollary 1 of the manuscript holds (expected
#' utility strictly decreasing in toxicity-efficacy correlation).
#'
#' @param alpha Efficacy exponent, \eqn{\alpha \ge 1}.
#' @param beta Toxicity exponent, \eqn{\beta \ge 1}.
#' @param lam Toxicity penalty rate, \eqn{\lambda > 0}.
#' @return A function \code{u(x_T, x_E)} returning \eqn{U(x_T, x_E)}.
#'   The returned function is vectorized: passing equal-length vectors
#'   for \code{x_T} and \code{x_E} returns a numeric vector of utilities.
#'
#' @details
#' \code{make_utility(1, 1, 2)} reproduces the utility used for
#' scenarios S1 to S6 of the manuscript; \code{make_utility(1, 2, 3)}
#' reproduces the utility used for the S7 block.
#'
#' @seealso \code{\link{utility_substitution}} for a single-call
#'   evaluation without constructing a closure.
#' @examples
#' u <- make_utility(1, 1, 2)
#' u(x_T = 0.2, x_E = 0.7)
#' u(x_T = c(0.1, 0.3), x_E = c(0.5, 0.8))
#' @export
make_utility <- function(alpha = 1, beta = 1, lam = 2) {
  stopifnot(alpha >= 1, beta >= 1, lam > 0)
  force(alpha); force(beta); force(lam)
  function(x_T, x_E) {
    x_E^alpha * exp(-lam * x_T^beta)
  }
}

#' Evaluate the substitution utility directly
#'
#' Convenience function equivalent to
#' \code{make_utility(alpha, beta, lam)(x_T, x_E)} without constructing
#' a closure.
#'
#' @param x_T Toxicity value(s) on the bounded scale \eqn{(0,1)}.
#' @param x_E Efficacy value(s) on the bounded scale \eqn{(0,1)}.
#' @inheritParams make_utility
#' @return Numeric vector of utility values.
#' @examples
#' utility_substitution(x_T = 0.1, x_E = 0.8, alpha = 1, beta = 1, lam = 2)
#' @export
utility_substitution <- function(x_T, x_E, alpha = 1, beta = 1, lam = 2) {
  x_E^alpha * exp(-lam * x_T^beta)
}
