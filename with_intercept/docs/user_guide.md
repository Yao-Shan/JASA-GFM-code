# Fitting an existing network

The model includes a distinct additive effect for every node. Supply a complete,
symmetric numeric matrix. Its diagonal is ignored. Poisson edges are nonnegative
integers; binary models require 0/1. Gaussian fits use a known positive `sigma`.
Missing edges, covariates and directed observations are not handled by this interface.

From the `with_intercept/` folder:

```r
source("R/gfm_baseline.R")
A <- as.matrix(read.csv("network.csv", header = FALSE))
fit <- gfb_fit(A, K = 2, model = "poisson", control = list(tolerance = 0.05))
print(fit)
str(fit$diagnostics)
saveRDS(fit, "fit.rds")
```

`gfb_fit` accepts model aliases `gaussian`, `logit`, `probit`, `poisson`.
Its output includes `U_hat`, `alpha_hat`, `Lambda_hat`, the observed matrix `A`,
`config`, and `diagnostics`. Fitting does not call an independent second optimizer
for inference. Inspect convergence, scaled-gradient residual, inequality margins,
and one-sided observed nodes. Failed or boundary fits should remain counted.

`U_hat` is centered and whitened. Its orientation is not identified without an
anchor. Simulation inference aligns to the known `U_star`; this known-truth step
must not be silently replaced by `U_hat` in a real-data claim.

## Inference from the same fit

When simulation truth is available:

```r
inference <- gfb_infer(fit, U_star = truth$U_star,
  alpha_star = truth$alpha_star, Lambda_star = truth$Lambda_star,
  nodes = 1, debias = TRUE, debias_method = "pdf_bp_bm")
```

`gfb_infer` reuses the supplied fit. It does not optimize the likelihood again.
Global debiasing requires the true `U_star` in this implementation. Without known
truth one may fit and predict; this package does not claim theorem-calibrated
coordinate errors or a validated truth-free global correction.
The T1 vector contains joint coordinates `u1`, `u2`, and `alpha_joint` at K=2;
its last coordinate is not the separately reported marginal-alpha statistic.
T2 has raw and corrected symmetric-Lambda coordinates. `correction_score_ell2`
and `correction_BM` are aliases for one term, not two corrections to add.

To defer inference, retain `fit.rds` and the matching truth/configuration. Later
read that fit and call `gfb_infer`; do not generate a new fit or substitute a new
network. Production studies additionally verified network/fit identity hashes.

## Computational controls

Core defaults are intentionally preserved. In particular the low-level
`gfb_config` tolerance is 0.001; completed G1-G5 studies explicitly used 0.05.
The study configuration helper specifies the recorded controls rather than silently
changing library defaults. The resource allowance is 2048 MiB for n=1000 and
6144 MiB for n=2000 in those studies. It is an admission bound, not a promise that
that much RAM suffices for many concurrent workers.

This upload package uses ordinary R for fitting and inference. The simulation
driver defaults to this implementation, and the notebooks do not load a runtime
folder or require Rcpp. For parallel simulation, use the standard PSOCK backend
and an explicit worker count as described in the simulation guide.
