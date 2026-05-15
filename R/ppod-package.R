#' ppod: Posterior Predictive Optimal Biological Dose Selection
#'
#' A Bayesian phase I/II adaptive design for optimal biological dose (OBD)
#' selection with continuous toxicity and efficacy endpoints. Dose
#' assignment and final selection are based on posterior predictive
#' expected utility under the joint toxicity-efficacy distribution on a
#' bounded transformed scale.
#'
#' @section Main entry points:
#' \describe{
#'   \item{\code{\link{simulate_trial}}}{Simulate one adaptive trial under a
#'     given scenario and decision criterion (joint, marginal, or plug-in).}
#'   \item{\code{\link{run_replicates}}}{Run many replicates of a
#'     (scenario, design) cell and collect dose selections.}
#'   \item{\code{\link{operating_characteristics}}}{Compute correct OBD
#'     selection probability, mean utility loss, near-optimal selection
#'     rate, and patient allocation from a replicate set.}
#'   \item{\code{\link{boin12_trial}}}{Binary-endpoint comparator
#'     simulator (Lin et al. 2020) used for benchmarking.}
#'   \item{\code{\link{true_expected_utility}}}{Deterministic-truth
#'     expected utility per dose under a scenario.}
#'   \item{\code{\link{ppod_scenarios}}}{The 12 simulation scenarios from
#'     the manuscript, available as a built-in dataset.}
#' }
#'
#' @section Decision criteria:
#' \describe{
#'   \item{\code{\link{eu_joint}}}{Posterior predictive expected utility
#'     under the joint Gaussian working model.}
#'   \item{\code{\link{eu_marginal}}}{Posterior predictive expected utility
#'     evaluated as if toxicity and efficacy were independent within dose.}
#'   \item{\code{\link{eu_plugin}}}{Utility at posterior mean dose-response
#'     estimates; ignores residual variability and parameter uncertainty.}
#' }
#'
#' @section Inference:
#' \describe{
#'   \item{\code{\link{fit_laplace}}}{Laplace approximation to the joint
#'     posterior over softplus toxicity slopes, random-walk efficacy
#'     deviations, log-scale residual standard deviations, and atanh
#'     correlation.}
#'   \item{\code{\link{sample_posterior}}}{Multivariate normal sampler from
#'     the Laplace posterior.}
#' }
#'
#' @docType package
#' @name ppod-package
#' @aliases ppod
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
