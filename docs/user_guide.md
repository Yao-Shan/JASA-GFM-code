# Single-network estimation guide

Use `R/gfm_estimation.R` to fit one network. It defines functions only and
requires base R; it does not simulate data or run Monte Carlo repetitions.

**Lambda bias correction requires the true latent-coordinate matrix `U_star`.**
When truth is unknown, use `debias = FALSE` to obtain ordinary estimates.
The function does not substitute estimated coordinates for truth.

## 1. Prepare the inputs

| Input | Requirement |
| --- | --- |
| `A` | Numeric, symmetric `n × n` network matrix with every off-diagonal entry observed. The diagonal is ignored. |
| `K` | Latent dimension, an integer with `1 ≤ K < n`. All output dimensions adapt to it. |
| `model` | `"logit"`, `"probit"`, `"poisson"`, or `"gaussian"`. |
| `U_star` | Required for correction: the true `n × K` matrix for **all nodes**, in the same node order as `A`, satisfying `crossprod(U_star) / n = diag(K)`. |
| `sigma` | Known positive noise standard deviation, required only for Gaussian observations. |

Logit and Probit require binary edges; Poisson requires nonnegative integer
counts; Gaussian accepts real-valued edges. Missing edges, directed networks,
covariates, and a separate intercept are not supported.

Run commands from the repository root. For CSV files without headers or row labels:

```r
source("R/gfm_estimation.R")
A <- gfm_read_matrix("network.csv")
U_star <- as.matrix(read.csv("U_star.csv", header = FALSE))
```

For files with column headers and node labels in the first column:

```r
A <- gfm_read_matrix("network_named.csv", node_names = TRUE)
U_star <- as.matrix(read.csv("U_star_named.csv", row.names = 1,
                            check.names = FALSE))
```

When both inputs have node labels, their order must match exactly. Without
labels, rows are matched by position. The function checks dimensions, finite
values, and normalization (maximum absolute discrepancy `1e-6`); it does not
reorder, center, or normalize `U_star`. Supply the true coordinates used by
your model, not `U_hat` or an arbitrary normalized matrix. These numerical
checks cannot establish whether supplied values are genuinely the truth.

## 2. Fit the model

With known true coordinates:

```r
fit <- gfm_fit(A, K = ncol(U_star), model = "logit", U_star = U_star)
print(fit)
fit$diagnostics
fit$debiasing
```

Change `model` for the other observation types. For Gaussian data, also supply
its known noise standard deviation, for example:

```r
fit <- gfm_fit(A, K = ncol(U_star), model = "gaussian",
               sigma = 0.5, U_star = U_star)
```

If `U_star` is unknown:

```r
fit <- gfm_fit(A, K = 2, model = "logit", debias = FALSE)
```

This returns ordinary estimates only. Requesting correction without `U_star`
raises an error before optimization. The true `Lambda_star` is not needed
for correction; it is needed only when evaluating estimation errors against
truth in a simulation.

## 3. Read and save the results

| Output | Meaning |
| --- | --- |
| `fit$U_hat` | Estimated node coordinates, `n × K`. |
| `fit$Lambda_hat` | Ordinary symmetric parameter estimate, `K × K`. |
| `fit$Lambda_debiased` | Bias-corrected Lambda; `NULL` unless correction succeeds. |
| `fit$diagnostics` | Optimizer convergence, iteration count, and numerical checks. |
| `fit$debiasing` | Correction status, explanation, and coordinate frame. |
| `fit$alignment` | Rotation to the supplied truth and alignment diagnostics; `not_requested` when truth is omitted. |

These are parameter estimates, not standardized simulation errors. The
`lambda_v_*` fields store diagonal entries divided by `sqrt(2)`, followed by
upper-triangle off-diagonal entries in pair order `(1,2), (1,3), …`.

Correction status is `success`, `failed`, `skipped` (fit did not converge), or
`not_requested`. A correction failure retains the ordinary estimates and emits
a warning. Inspect the status before using a corrected result.

```r
saveRDS(fit, "gfm_fit.rds")
if (identical(fit$debiasing$status, "success")) {
  write.csv(fit$Lambda_debiased, "Lambda_debiased.csv")
  fitted_means <- gfm_predict(fit, debiased = TRUE)
}
```

`gfm_predict(fit)` uses ordinary Lambda; `type = "link"` returns the fitted
linear predictor. Predictions require a converged fit, preserve node order,
and have `NA` on the diagonal. No confidence intervals are returned.

## 4. Adjust the controls

The model-specific optimizer defaults are:

| Model | `max_iter` | `learning_rate` | `tolerance` |
| --- | ---: | ---: | ---: |
| Logit | 40000 | 1 | 0.001 |
| Probit | 40000 | 1 | 0.001 |
| Gaussian | 10000 | 0.1 | 0.001 |
| Poisson | 60000 | 0.02 | 0.001 |

For example, to allow more iterations and tighten the stopping tolerance:

```r
fit <- gfm_fit(A, K = ncol(U_star), model = "logit", U_star = U_star,
               control = list(max_iter = 60000L, tolerance = 0.01))
```

`max_iter` counts optimization updates for this one network, not simulated
networks. If convergence fails, inspect the diagnostics before increasing it;
a smaller `learning_rate` may help unstable updates. A smaller `tolerance`
requests a stricter stopping criterion and can take longer.

Other defaults are `M_bound_U = 20`, `M_bound_Lambda = 20`,
`probability_eps = 0.001`, and `max_debias_memory_mb = 2048`. The first two
set parameter constraints (`M_bound_U ≥ sqrt(K)`); `probability_eps` controls
Bernoulli probability clipping. Change them only when appropriate for the model.

Correction uses dense matrices of dimension `n*K`. Its memory guard uses a
planning estimate of `16 * 8 * (n*K)^2 / 1024^2` MiB, not a peak-memory guarantee.
Raise the limit only if sufficient memory is available. The internal
pseudoinverse helper removes exactly `n_drop` eigen-directions and requires all
retained eigenvalues to be positive; it has no adjustable cutoff tolerance.

## 5. How truth enters the correction

Optimization uses `A` and the model settings. If `U_star` is supplied, the
result is then aligned to its coordinates using the same Procrustes rotation
as the simulation notebook:

$$
R = \arg\min_{Q^\top Q=I_K}\|\widehat U Q-U^\star\|_F,
\qquad
\widehat U_{\mathrm{a}}=\widehat U R,
\qquad
\widehat\Lambda_{\mathrm{a}}=R^\top\widehat\Lambda R.
$$

Both returned ordinary estimates use this aligned frame, even when
`debias = FALSE` if truth was supplied. Alignment preserves the fitted signal
`U_hat %*% Lambda_hat %*% t(U_hat)`. Without supplied truth, ordinary estimates
remain in the optimizer's coordinate frame.

The correction evaluates the information and bias terms at the aligned
estimates, with **`U_anchor = U_star`** in the constraint derivatives. It adds
the resulting correction to aligned Lambda and does not bias-correct U.
This follows the simulation's use of true coordinates; it is not a correction
available from an observed network alone when those coordinates are unknown.

To verify the implementation locally, run `Rscript --vanilla tests/test_gfm.R`.
The checks compare alignment and correction with the simulation formulas,
cover all four likelihoods and adaptive K, and do not run Monte Carlo experiments.
