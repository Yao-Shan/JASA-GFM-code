# JASA GFM code

R implementations of a symmetric generalized factor model for Gaussian,
Poisson, Bernoulli-logit, and Bernoulli-probit network observations.

The repository provides an interface for fitting **one existing network** and
a separate simulation notebook used as the algorithm reference.

## Analyze an existing network

Run the following from the repository root:

```r
source("R/gfm_estimation.R")

# Your observed square network matrix, without a header or row labels.
A <- gfm_read_matrix("network.csv")
fit <- gfm_fit(A, K = 2, model = "logit")

fit$diagnostics        # Convergence and numerical diagnostics
fit$debiasing          # Whether the Lambda correction succeeded
fit$U_hat             # Estimated coordinates of all observed nodes
fit$Lambda_hat        # Raw symmetric parameter estimate
fit$Lambda_debiased   # Corrected estimate; NULL when unavailable

if (identical(fit$debiasing$status, "success")) {
  fitted_means <- gfm_predict(fit, debiased = TRUE)
}
```

Use `model = "probit"`, `"poisson"`, or `"gaussian"` for the other models.
Gaussian fits require an explicitly supplied, known noise standard deviation,
for example `sigma = 1`. All parameter dimensions adapt to the chosen `K`;
the function does not select `K` automatically.

The independent estimation source requires only base R. It does not generate
data, set a random seed, start Monte Carlo repetitions, or write output files.
Optimization iterations are updates for the single supplied dataset.

Read the [practitioner guide](docs/user_guide.md) for named-node CSVs,
controls, returned values, failures, prediction, and saving results.

## Repository structure

```text
JASA-GFM-code/
├── R/
│   └── gfm_estimation.R       # Functions for one observed network
├── docs/
│   └── user_guide.md          # Inputs, controls, results, and interpretation
├── simulations/
│   └── node_inference.Rmd     # Data generation and Monte Carlo analysis
├── tests/
│   └── test_gfm.R             # Deterministic regression checks
└── README.md
```

For an existing dataset, start with [`R/gfm_estimation.R`](R/gfm_estimation.R)
and the [user guide](docs/user_guide.md). To run the simulation study, use
[`simulations/node_inference.Rmd`](simulations/node_inference.Rmd).

## Input and interpretation

Supply both triangles of a symmetric numeric matrix with every off-diagonal
edge observed. Binary models require 0/1 entries; Poisson requires nonnegative
integer counts. The diagonal is ignored. The model has no separate intercept,
covariate, missing-edge, or directed-network interface.

Latent coordinates are identifiable only up to orthogonal transformations.
The practical correction uses `U_anchor = U_hat` and returns Lambda in the
fitted coordinate frame. Simulation diagnostics can instead align to known
truth. The practical interface does not report theorem-calibrated coordinate
confidence intervals for its data-dependent frame, and it does not debias U.

Always inspect convergence and correction status. The implementation uses
dense matrices, including information matrices of dimension `n*K` by `n*K`.
The usage guide describes the memory guard and computational limits.

## Verification

```sh
Rscript --vanilla tests/test_gfm.R
```

The current suite has 80 checks covering all four likelihoods, agreement with
the simulation formulas under the same coordinate anchor, adaptive dimensions,
rotation invariance, input validation, failure handling, and predictions.
Fixtures are used only inside the test script. The suite extracts selected
function definitions from the notebook without running its simulation pipeline.

## Simulation reference

Open `simulations/node_inference.Rmd` in RStudio to reproduce or adapt the simulation
workflow. It uses `rmarkdown`, `knitr`, `doParallel`, and `foreach`, which are
not required by the practitioner interface or the test suite.

Review the YAML parameters before knitting: `run_simulation: true` starts the
configured experiment; the current notebook specifies 1,000 replications and
up to seven workers. Set `run_simulation: false` for its truth-checking mode.
The notebook writes its configured results under `data/` at the repository
root. Simulation results and analysis datasets are generated or supplied
separately; they are not bundled in this code repository.
