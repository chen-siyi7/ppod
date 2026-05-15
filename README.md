# ppod: Posterior Predictive Optimal Biological Dose Selection

<!-- badges: start -->
[![R-CMD-check](https://github.com/author/ppod/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/author/ppod/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

`ppod` is an R package implementing a Bayesian phase I/II adaptive
design for optimal biological dose (OBD) selection with **continuous**
toxicity and efficacy endpoints. Dose assignment and final selection
are based on posterior predictive expected utility under the joint
toxicity-efficacy distribution on a bounded transformed scale.

This is the companion software for the manuscript:

> *Continuous-Endpoint Posterior Predictive Utility with
> Dependence-Sensitive Benefit-Risk Evaluation for Optimal Biological
> Dose Selection* (under review, 2026).

## What the package provides

- Three posterior predictive expected-utility criteria for OBD selection:
  joint (`eu_joint`), marginal (`eu_marginal`), and plug-in
  (`eu_plugin`).
- A monotone toxicity dose-response model and a second-order
  random-walk efficacy model that admits monotone, plateau, and
  unimodal patterns.
- Posterior predictive admissibility gates aligned with the utility
  criterion (`admissibility`).
- Adaptive trial simulation under the full design (`simulate_trial`),
  including the run-in escalation rule and the no-skipping safeguard.
- A simplified BOIN12 binary-endpoint comparator (`boin12_trial`)
  matching Lin et al. (2020).
- The 12 simulation scenarios from the manuscript as built-in data
  (`ppod_scenarios`), including the S4d dose-varying dependence stress
  test and the S7 close-tied family.
- Operating-characteristics summarization
  (`run_replicates`, `operating_characteristics`).
- Reproducibility scripts in `inst/scripts/` that reproduce the
  main-text simulation tables and the supplementary sensitivity
  studies.

## Installation

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("author/ppod")
```

The package has no compiled code and depends only on `MASS` and base R.

## Quick start

```r
library(ppod)

# Inspect a scenario
scen <- ppod_scenarios$S4b
str(scen)

# Run a single adaptive trial under the joint posterior predictive
# utility criterion.
set.seed(20260515)
result <- simulate_trial(scen, criterion = "joint")
result$selected     # selected OBD index
table(result$dose_history)

# Compare three continuous-endpoint criteria on a small replicate set.
designs <- c("joint", "marginal", "plugin")
cells <- lapply(designs, function(d) {
  run_replicates(scen, design = d, n_reps = 20L, seed_base = 1000L)
})
names(cells) <- designs
sapply(cells, function(c) operating_characteristics(c, scen)$p_correct)
```

## Reproducing the manuscript

The scripts in `inst/scripts/` reproduce the manuscript figures and
tables.

```r
# Path to the bundled scripts
sys.file <- system.file("scripts", package = "ppod")
list.files(sys.file)
# reproduce_main.R reproduces Tables 3 and Figures 1, 2 of the main text
# reproduce_supp.R reproduces the four supplementary sensitivity studies
```

These are full simulation studies: `reproduce_main.R` runs 400
replicates across 12 scenarios and 4 designs (approximately 30 to 60
minutes on a recent laptop), and `reproduce_supp.R` runs four
secondary studies (approximately 60 to 90 minutes total).

## Decision criteria at a glance

For a trial with accumulated data \(\mathcal{D}_n\), the package
implements three utility criteria at dose \(d_j\):

| Criterion  | Formula                                                                | Captured uncertainty |
|------------|------------------------------------------------------------------------|----------------------|
| `eu_joint` | \(E[U(g_T(Y_T), g_E(Y_E)) \mid d_j, \mathcal{D}_n]\) under joint model | parameter + residual + dependence |
| `eu_marginal` | Same integral with \(\rho_j \equiv 0\) at the utility step          | parameter + residual (no dependence) |
| `eu_plugin` | \(U(g_T(\hat\mu_T(d_j)), g_E(\hat\mu_E(d_j)))\) at posterior mode    | none |

The joint-versus-marginal contrast isolates the role of
dependence-sensitive integration; the joint-versus-plug-in contrast
isolates posterior predictive integration over parameter and residual
uncertainty.

## Posterior inference

Posterior inference uses a Laplace approximation
(`fit_laplace`, `sample_posterior`) suitable for large
operating-characteristics studies where full MCMC at every cohort would
be prohibitive. A Laplace-versus-MCMC sensitivity analysis in the
supplement of the manuscript shows that Laplace is accurate in
well-separated scenarios but understates the joint design's ability to
recover dose-rank differences in close-tied scenarios with informative
dependence; users running a single trial in practice should consider
full MCMC posterior sampling and pass the resulting samples directly to
`eu_joint`.

## Citing

If you use `ppod` in published work, please cite the accompanying
manuscript:

```bibtex
@article{author2026ppod,
  title   = {Continuous-Endpoint Posterior Predictive Utility with
             Dependence-Sensitive Benefit-Risk Evaluation for Optimal
             Biological Dose Selection},
  author  = {Author One and Author Two and Author Three},
  journal = {Biometrics},
  year    = {2026},
  note    = {Under review}
}
```

## License

MIT (see `LICENSE.md`).

## Contributing

Issues and pull requests are welcome at
<https://github.com/author/ppod>. Before submitting changes, please
run:

```r
devtools::document()
devtools::check()
devtools::test()
```

## Acknowledgments

We thank colleagues who provided feedback on the methodology and
implementation.
