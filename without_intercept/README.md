# GFM without intercepts

Fit symmetric generalized factor models for Gaussian, Poisson, Logit, and
Probit networks. Run the commands below from this `without_intercept/` folder.
For node intercepts, use the [separate implementation](../with_intercept/README.md).

## Choose a workflow

| Task | Entry point | Instructions |
| --- | --- | --- |
| Fit one observed network | [R/gfm_estimation.R](R/gfm_estimation.R) | [Estimation guide](docs/user_guide.md) |
| Run Monte Carlo simulations | [simulations/node_inference.Rmd](simulations/node_inference.Rmd) | [Simulation guide](docs/simulation_guide.md) |
| Draw summary QQ plots, sample-size comparisons, or all-node figures | [plots/plot_figures.Rmd](plots/plot_figures.Rmd) | [Simulation guide](docs/simulation_guide.md) |
| Export individual QQ plots and Lambda histograms | [plots/plot_individual.Rmd](plots/plot_individual.Rmd) | [Individual figures](docs/simulation_guide.md#individual-qq-plots-and-lambda-histograms) |

## Fit an observed network

Supply a symmetric numeric matrix with all off-diagonal entries observed.
For a CSV without headers or row labels:

```r
source("R/gfm_estimation.R")
A <- gfm_read_matrix("network.csv")
fit <- gfm_fit(A, K = 2, model = "logit", debias = FALSE)
fit$diagnostics
fit$U_hat
fit$Lambda_hat
```

**Lambda bias correction requires the true coordinates `U_star` for all nodes.**
When they are unknown, use `debias = FALSE` as above. This returns ordinary
estimates; it does not substitute estimated coordinates for truth. See the
[estimation guide](docs/user_guide.md) for correction, input requirements,
optimizer controls, and returned values. Gaussian fits also require a known
noise standard deviation, supplied through `sigma`.

## Simulations and saved figures

The simulation notebook defaults to `run_simulation: true`, `n: 1000`, and
`n_sim: 1000`. To inspect its settings and construct the truth without running
network replications:

```r
rmarkdown::render("simulations/node_inference.Rmd",
                  params = list(run_simulation = FALSE))
```

The [simulation guide](docs/simulation_guide.md) explains paper settings,
small trial runs, parameter changes, output paths, and plotting requirements.
Both plotting notebooks start with `run_plots: false` and read saved results
without fitting models. Summary figures and individual figures have separate
entries because they serve different layouts. Theorem 2 QQ plots show only
bias-corrected statistics.

Saved outputs and reports are retained locally pending review. To draw QQ plots
and histograms, supply matching replication-level files or generate them with
the simulation notebook.

## Checks and provenance

```sh
Rscript --vanilla tests/test_gfm.R
Rscript --vanilla tests/test_simulation_config.R
Rscript --vanilla tests/test_plotting.R
```

These check estimation formulas and interfaces, simulation configuration, and
plotting data handling; they do not run the paper's Monte Carlo study.
See [source provenance](docs/provenance.md) for the code snapshot and the
relationship to historical results.
