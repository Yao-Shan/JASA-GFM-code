# GFM with additive node effects

This folder implements `Theta_ij = alpha_i + alpha_j + u_i^T Lambda u_j` for Gaussian,
Bernoulli-logit, Bernoulli-probit and Poisson networks. It is organized in the same
way as the sibling `without_intercept/` folder. Fitting and inference formulas
are preserved from the completed studies. The simulation driver defaults to
ordinary R; local acceleration sources are excluded from this upload package.

## Start here

Run the following with this folder as your working directory:

```r
source("R/gfm_baseline.R")
A <- as.matrix(read.csv("network.csv", header = FALSE))
fit <- gfb_fit(A, K = 2, model = "poisson", control = list(tolerance = 0.05))
fit$diagnostics
fit$alpha_hat
fit$U_hat
fit$Lambda_hat
# Only for a converged fit:
# fitted_means <- gfb_predict(fit)
```

For an observed network, read [the user guide](docs/user_guide.md). The global
simulation correction requires the true `U_star`; it is not a truth-free
practitioner correction. For simulation, open
[`simulations/node_inference.Rmd`](simulations/node_inference.Rmd). Its default
`run_simulation: false` inspects a small design without fitting networks.
The all-node entry is [`simulations/all_nodes.Rmd`](simulations/all_nodes.Rmd).

| Folder | Contents |
|---|---|
| `R/` | Four numerical modules; source `gfm_baseline.R` |
| `simulations/` | Small trial notebooks and explicit G1-G5 study settings |
| `plots/` | Functions and notebooks for QQ plots, Lambda histograms and summary plots |
| `docs/` | Inputs, seeds, study protocol, implementation limitations |
| `tests/` | Package loading, settings, formula and small execution checks |

See [the simulation guide](docs/simulation_guide.md), [seed definitions](docs/seeds.md),
and [implementation notes](docs/implementation_notes.md).
Saved outputs and reports are retained locally pending review. To create figures,
supply matching saved records or generate them with the simulation entry points.

## Check the portable package

```sh
Rscript --vanilla tests/test_package.R
Rscript --vanilla tests/test_small_run.R
Rscript --vanilla tests/test_runtime_free.R
Rscript --vanilla tests/test_plots.R
```

The first check does not fit networks. The third compares ordinary-R sequential
and PSOCK execution without acceleration sources. The second runs four small deterministic-seed
execution checks; it is not a Monte Carlo calibration study. Neither command accesses
local historical directories. R base/stats/utils suffice for numerical use. Notebooks
require `knitr` and `rmarkdown`; Pandoc is needed to render HTML. Rcpp and a C++
compiler are not required by these notebooks. Standard PSOCK parallel execution
remains available. New simulation output is generated locally and excluded from
the upload package.

For figures, open [`plots/plot_figures.Rmd`](plots/plot_figures.Rmd). This single
entry reads saved results and exports QQ plots or either Lambda histogram type.
Plotting is disabled by default; see the [simulation guide](docs/simulation_guide.md#plot-saved-results)
for input paths and plot selection. Theorem 2 QQ plots show debiased statistics only.
