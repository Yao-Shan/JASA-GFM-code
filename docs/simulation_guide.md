# Simulation and plotting guide

Use [`simulations/node_inference.Rmd`](../simulations/node_inference.Rmd) to
generate networks, fit the model repeatedly, and assess the nodewise and Lambda
statistics. Use [`plots/plot_figures.Rmd`](../plots/plot_figures.Rmd) to draw figures
from results that have already been saved. The plotting notebook does not fit
models or run simulations.

For one observed network, use the separate [practitioner guide](user_guide.md).

## Setup and path conventions

Run the R commands in this guide from the repository root. The simulation
notebook requires `rmarkdown`, `knitr`, `doParallel`, and `foreach`:

```r
install.packages(c("rmarkdown", "knitr", "doParallel", "foreach"))
```

Install these packages only if they are missing. Knitting HTML also requires
Pandoc, which is included with RStudio. Plot calculations use base R;
`rmarkdown` and `knitr` provide its HTML notebook.

During rendering, paths in notebook parameters are relative to the notebook's
directory. For example, `../data/node1` means the repository's `data/node1`
directory from either `simulations/` or `plots/`. Absolute paths are also
accepted.

## Paper simulation settings

The simulation notebook's defaults run the main setting below. Knitting it
starts the full experiment because `run_simulation` is `true`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `model` | `"all"` | Gaussian, Poisson, Logit, and Probit |
| `n` | `1000` | Number of nodes in each simulated network |
| `n_sim` | `1000` | Number of independent network replications conditional on the latent design |
| `lambda_values` | `"0.5,0.5"` | Diagonal entries of the true Lambda |
| `k` | `NULL` | Infer K from the number of supplied Lambda entries |
| `target_node` | `1` | Node used for Theorem 1 diagnostics |
| `design_seed` | `238687` | Seed for constructing the latent design |
| `replication_seed_base` | `2000` | Replication r uses seed `2000 + r` |
| `n_cores` | `7` | Requested maximum number of parallel workers |
| `truth_check_action` | `"warning"` | Report regularity-threshold issues without stopping on those issues |
| `output_dir` | `NULL` | Derive a destination from the actual simulation configuration |

With `model = "all"`, models run sequentially; replications within a model use
parallel workers. The notebook caps the worker count at seven and may use fewer
workers when appropriate. The Gaussian simulation has known noise standard
deviation `sigma = 1`. Each model retains its original learning rate, maximum
iteration count, convergence tolerance, bounds, and numerical thresholds when
the corresponding YAML overrides are `null`.

To run the default main setting:

```r
rmarkdown::render("simulations/node_inference.Rmd")
```

To inspect the constructed truth and its diagnostics without running network
replications:

```r
rmarkdown::render(
  "simulations/node_inference.Rmd",
  params = list(run_simulation = FALSE),
  output_file = "truth_check.html"
)
```

This mode constructs and checks the latent truth. It does not load an earlier
result CSV or redraw its figures.

### Run both Lambda settings and all three sample sizes

The paper settings use diagonal Lambda values `(0.5, 0.5)` and `(0.7, -0.5)`,
with `n = 250`, `500`, and `1000`. The following explicitly runs all six settings,
with all four models and 1,000 replications per model in each setting:

```r
lambda_settings <- c(equal = "0.5,0.5", mixed = "0.7,-0.5")
for (n_value in c(250L, 500L, 1000L)) {
  for (setting_name in names(lambda_settings)) {
    rmarkdown::render(
      "simulations/node_inference.Rmd",
      params = list(
        model = "all",
        run_simulation = TRUE,
        n = n_value,
        n_sim = 1000L,
        k = NULL,
        lambda_values = unname(lambda_settings[setting_name]),
        target_node = 1L,
        design_seed = 238687L,
        replication_seed_base = 2000L,
        n_cores = 7L,
        output_dir = NULL
      ),
      output_file = paste0("node_inference_n", n_value, "_", setting_name, ".html")
    )
  }
}
```

This is a substantial Monte Carlo experiment. To run just one setting, use one
`render()` call with its chosen `n` and `lambda_values`.

The latent design is fixed across replications within a setting. The design seed
`238687` was selected after inspecting simulation performance; it is not an
independently selected validation seed. Changing `replication_seed_base` changes
the generated networks while retaining the same latent design. Changing `n`
can change the realized latent coordinates even with the same design seed;
node 1 remains the same row index, but is not guaranteed the same coordinate
vector across sample sizes.

These commands reproduce the configured experiment. The repository does not
bundle the saved paper result CSVs. For exact redrawing of an existing figure,
use the original saved results and metadata. Numerical results may also depend
on the R and numerical-library environment; execution settings are saved with
each run.

## A small trial run

Reduce network size and replication count when checking that the workflow runs
on your computer. For example:

```r
rmarkdown::render(
  "simulations/node_inference.Rmd",
  params = list(
    model = "gaussian",
    run_simulation = TRUE,
    n = 50L,
    n_sim = 5L,
    k = NULL,
    lambda_values = "0.5,0.5",
    target_node = 1L,
    n_cores = 1L,
    truth_check_action = "warning",
    output_dir = "../data/trial/gaussian"
  ),
  output_file = "trial_gaussian.html"
)
```

This example retains the default seeds and model-specific optimizer settings.
Change `model` to another model, or to `"all"` with an output directory such as
`"../data/trial"`, to try the other likelihoods. A small trial checks execution
and output structure; five replications do not establish normal approximation,
coverage, or other statistical performance. Small networks can produce truth
diagnostic warnings or estimation failures. Inspect these diagnostics before
interpreting a trial result.

The main controls have different roles:

| Control | What changes |
| --- | --- |
| `n` | The size of each simulated network and the latent design |
| `n_sim` | How many networks are generated and fitted; fewer repetitions make Monte Carlo summaries less precise |
| `max_iter` | The upper limit on optimizer updates for each single network; reducing it can cause nonconvergence |
| `n_cores` | Parallel execution, rather than the intended statistical setting |
| `lambda_values` | The true Lambda and, when `k = NULL`, the latent dimension |
| `target_node` | Which node's Theorem 1 statistics are saved |

Start by reducing `n_sim` and `n`, and by choosing a suitable worker count.
Do not reduce `max_iter` merely to imitate a smaller dataset. Keep the original
learning rates, tolerances, and numerical thresholds for paper comparisons.

### Change the latent dimension

For a three-dimensional simulation, supply three diagonal entries:

```r
params = list(k = NULL, lambda_values = "0.5,0.4,-0.3")
```

This is an example parameter fragment for a `render()` call. The notebook then
uses K = 3. If `k` is explicitly supplied, it must equal the number of Lambda
entries. Both the nodewise K components and the K(K+1)/2 components of a
symmetric Lambda adapt to this dimension. This is dimension-aware computation,
not a procedure for choosing K from data.

## Saved results

With `output_dir = NULL`, the notebook derives a directory of the form

```text
data/node{node}/n{n}_k{k}_lambda_{value-tag}_target{node}/{model}/
```

The Lambda tag contains readable values and a short digest. File names also
record replication counts, seeds, and numerical settings. Re-running exactly
the same configuration in the same destination replaces its matching output
files. Use a separate explicit output directory when retaining another run.

For an explicit `output_dir`, one selected model writes directly into that
directory. With `model = "all"`, a model-named subdirectory is appended.

| File ending | Contents |
| --- | --- |
| `_results.csv` | One row per attempted replication, standardized statistics, and diagnostics, including failed rows |
| `_config.rds` | Model, sample size, Lambda, seeds, and algorithm configuration |
| `_truth.rds` | The simulated latent truth |
| `_settings.rds` | Execution settings |
| `_complete.rds` | Results retaining their R attributes |
| `_qqplots.pdf` | Basic notebook QQ plots, when successful replications are available |

Keep result CSVs together with their matching configuration files. The names
`raw_lambda_*` and `debiased_lambda_*` in the simulation CSV refer to
**standardized errors**, rather than raw or corrected parameter estimates.
The plotting notebook uses these statistics as saved; it does not recenter or
rescale the observations to improve agreement with their reference distribution.

## Redraw figures without running simulations

Render the plotting notebook after simulation output is available:

```r
rmarkdown::render(
  "plots/plot_figures.Rmd",
  params = list(
    data_dir = "../data/node1",
    output_dir = "../data/figures"
  )
)
```

It searches the input directory recursively for native `_results.csv` files and
requires a matching `_config.rds` for each file. One setting is sufficient:
the notebook draws figures from the available combinations, without requiring
all models, both Lambda values, or all three sample sizes. Its optional
`sample_sizes` parameter selects sample sizes, and `main_n` selects the main QQ
sample size; by default each K/Lambda group uses its own largest available
sample size. If an explicit `main_n` is missing from a group, the notebook
reports that its main QQ figures were skipped. `max_columns` (default 4)
controls how many panels are placed across a page.

For example, redraw one small trial independently of the paper results:

```r
rmarkdown::render(
  "plots/plot_figures.Rmd",
  params = list(
    data_dir = "../data/trial/gaussian",
    output_dir = "../data/trial_figures",
    max_columns = 3L
  ),
  output_file = "trial_figures.html"
)
```

The plots use open circles and a dashed identity reference line, without a
reference band. Node components remain separate. Outputs include PDF and PNG
figures, the rendered HTML notebook, and `plot_inputs.csv`, `row_selection.csv`,
and `figure_manifest.csv` records.

For comparable panels, the notebook uses the same retained replications across
the node, raw Lambda, corrected Lambda, and joint statistics: a retained row
must have `convergence_code == 0` and finite values for all required statistics.
Attempted, converged, used, and excluded counts are reported, with the selected
original row indices. These plots describe the retained replications; failure
rates remain part of the experiment's results.

The reader rejects duplicate model/sample-size/Lambda cases and incompatible
seeds, target nodes, or truth-design settings. Point `data_dir` to a narrower
collection when comparing another experiment, rather than mixing independent
runs in one figure.

### Optional all-node figures

The nodewise simulation notebook evaluates the selected target node. It does
not produce the summaries for every node required by the all-node boxplots
and heatmaps.

If you already have native all-node output, pass its `component_summary.csv`
through `all_nodes_file`. The matching `config.rds` must be present in the same
directory. A directory or a vector of paths can also be supplied. For example:

```r
rmarkdown::render(
  "plots/plot_figures.Rmd",
  params = list(
    data_dir = "../data/node1",
    output_dir = "../data/figures",
    all_nodes_file = "../data/all_nodes/run_name/component_summary.csv"
  )
)
```

Replace `run_name` with an actual native all-node output directory. Each case
must contain complete summaries for every node and every component, with finite
summary values and replication counts matching the saved configuration.
Incomplete or mismatched summaries produce an error. Directory searches use
the per-case summaries and skip a run's aggregate copy when matching per-case
files are present. A legacy merged summary CSV without its corresponding
configuration metadata is insufficient.
Leave `all_nodes_file = NULL` to skip these figures.

## An existing smaller dataset

An observed dataset does not need a simulation setting. For example:

```r
source("R/gfm_estimation.R")
A <- gfm_read_matrix("network.csv")
fit <- gfm_fit(A, K = 2, model = "logit")
fit$diagnostics
fit$debiasing
fit$Lambda_debiased
```

Here `n` is inferred from the dimensions of `A`. There is no `n_sim`, design
seed, replication seed, or latent-truth input. The optimizer still performs
iterations to fit this one network. Consult the [practitioner guide](user_guide.md)
for model assumptions, optimizer controls, and interpretation of the correction.
