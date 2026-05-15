# ppod 0.1.0

Initial release accompanying the manuscript "Continuous-Endpoint Posterior
Predictive Utility with Dependence-Sensitive Benefit-Risk Evaluation for
Optimal Biological Dose Selection."

## Features

- Joint, marginal, and plug-in posterior predictive expected utility
  criteria for OBD selection with continuous toxicity and efficacy
  endpoints.
- Posterior predictive admissibility gates.
- Laplace approximation to the joint posterior; full MCMC is not
  required for OBD decisions in the working model.
- Binary-endpoint BOIN12 comparator implementation.
- 12 simulation scenarios from the manuscript, including the S4d
  dose-varying dependence stress test and the S7 close-tied family.
- Reproducibility scripts in `inst/scripts/`.
