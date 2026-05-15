#' Laplace approximation to the joint posterior
#'
#' Fits the joint dose-response model of the manuscript by L-BFGS-B
#' optimization of the negative log posterior, then constructs the
#' Laplace approximation: a multivariate normal centered at the
#' posterior mode with covariance equal to the inverse of the
#' numerical Hessian at the mode. Used for fast posterior approximation
#' inside the adaptive operating-characteristics simulation, where full
#' MCMC at every cohort is computationally prohibitive.
#'
#' @param Y_T Numeric vector of toxicity observations on the working
#'   scale, length \eqn{n}.
#' @param Y_E Numeric vector of efficacy observations on the working
#'   scale, length \eqn{n}.
#' @param doses Integer vector of dose indices (1-based), length \eqn{n}.
#' @param J Integer, number of dose levels.
#' @param theta_init Optional warm-start parameter vector. If
#'   \code{NULL}, defaults computed from \code{J} are used.
#' @param maxit Maximum L-BFGS-B iterations (default 200).
#' @return A list with components:
#'   \describe{
#'     \item{\code{ok}}{Logical, \code{TRUE} if optimization converged.}
#'     \item{\code{theta}}{Posterior mode (vector of length \eqn{2J + 3}).}
#'     \item{\code{cov}}{Posterior covariance matrix at the mode.}
#'   }
#'   When \code{ok} is \code{FALSE}, \code{theta} contains the
#'   warm-start (or default) initialization and \code{cov} is a
#'   well-conditioned default; downstream code should typically retry
#'   with more data or different initialization.
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 30
#' doses <- rep(1:3, each = 10)
#' Y_T <- rnorm(n, mean = c(-1.5, -0.8, -0.2)[doses], sd = 0.5)
#' Y_E <- rnorm(n, mean = c(-0.5, 0.5, 0.7)[doses], sd = 0.5)
#' fit <- fit_laplace(Y_T, Y_E, doses, J = 3)
#' fit$ok
#' fit$theta
#' }
#' @seealso \code{\link{sample_posterior}} to draw from the Laplace
#'   posterior; \code{\link{neg_log_post}} for the objective function.
#' @export
fit_laplace <- function(Y_T, Y_E, doses, J, theta_init = NULL, maxit = 200L) {
  P <- .n_params(J)
  if (is.null(theta_init)) {
    theta_init <- .theta_default(J)
  }
  stopifnot(length(theta_init) == P)

  fit <- try(
    stats::optim(
      par     = theta_init,
      fn      = neg_log_post,
      Y_T     = Y_T,
      Y_E     = Y_E,
      doses   = doses,
      J       = J,
      method  = "L-BFGS-B",
      control = list(maxit = maxit, factr = 1e9)
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error") || fit$convergence != 0L) {
    return(list(ok = FALSE, theta = theta_init, cov = diag(0.1, P)))
  }

  theta_hat <- fit$par
  H <- .numerical_hessian(
    neg_log_post, theta_hat,
    list(Y_T = Y_T, Y_E = Y_E, doses = doses, J = J)
  )
  H <- 0.5 * (H + t(H))

  cov_mat <- try(solve(H + diag(1e-5, P)), silent = TRUE)
  if (inherits(cov_mat, "try-error")) {
    cov_mat <- MASS::ginv(H + diag(1e-3, P))
  }
  cov_mat <- 0.5 * (cov_mat + t(cov_mat))
  ev <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  if (any(ev <= 0)) {
    cov_mat <- cov_mat + diag(1e-4 - min(ev), P)
  }

  list(ok = TRUE, theta = theta_hat, cov = cov_mat)
}

#' Draw samples from a Laplace posterior approximation
#'
#' Generates \code{S} multivariate normal samples \eqn{\theta^{(s)} \sim
#' \mathcal{N}(\hat\theta, \widehat\Sigma)} from a Laplace approximation
#' produced by \code{\link{fit_laplace}}. Uses a Cholesky factorization
#' with an eigen-decomposition fallback for marginally non-positive
#' covariance matrices.
#'
#' @param theta_hat Posterior mode (vector of length \eqn{2J + 3}).
#' @param cov Posterior covariance matrix.
#' @param S Number of samples to draw.
#' @return An \code{S}-by-\eqn{2J + 3} numeric matrix.
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 30
#' doses <- rep(1:3, each = 10)
#' Y_T <- rnorm(n, mean = c(-1.5, -0.8, -0.2)[doses], sd = 0.5)
#' Y_E <- rnorm(n, mean = c(-0.5, 0.5, 0.7)[doses], sd = 0.5)
#' fit <- fit_laplace(Y_T, Y_E, doses, J = 3)
#' draws <- sample_posterior(fit$theta, fit$cov, S = 200)
#' dim(draws)
#' }
#' @export
sample_posterior <- function(theta_hat, cov, S) {
  P <- length(theta_hat)
  L <- try(chol(cov), silent = TRUE)
  Z <- matrix(stats::rnorm(S * P), nrow = S)
  if (inherits(L, "try-error")) {
    e <- eigen(cov, symmetric = TRUE)
    Lfb <- e$vectors %*% diag(sqrt(pmax(e$values, 1e-8)))
    return(matrix(theta_hat, S, P, byrow = TRUE) + Z %*% t(Lfb))
  }
  matrix(theta_hat, S, P, byrow = TRUE) + Z %*% L
}
