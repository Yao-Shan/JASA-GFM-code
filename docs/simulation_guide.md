# Simulation guide

This guide explains the paper's simulation settings, how to reproduce the
experiments and summarize their results, and how to explore other settings.

- [`simulations/node_inference.Rmd`](../simulations/node_inference.Rmd): generate networks, fit models, and evaluate the Theorem 1 and Theorem 2 statistics.
- [`plots/plot_figures.Rmd`](../plots/plot_figures.Rmd): draw figures from saved results.

For an existing observed network, see the [data analysis guide](user_guide.md).

## 1. Simulation settings

The paper considers Gaussian, Poisson, Logit, and Probit models with
250, 500, or 1000 nodes. The true Lambda matrix is diagonal, with
diagonal entries (0.5, 0.5) or (0.7, −0.5). Each setting uses
1000 network replications. Theorem 1 is evaluated at node 1;
Theorem 2 concerns the global parameter Lambda and does not
depend on the selected node.

The current YAML uses the main paper sample size and repetition count.
Reduce `n` and `n_sim` for a small trial:

| Parameter | Default | Meaning |
| --- | --- | --- |
| `model` | `"all"` | Run all four models |
| `n` | `1000` | Number of nodes per network |
| `n_sim` | `1000` | Number of network replications |
| `lambda_values` | `"0.5,0.5"` | Diagonal entries of the true Lambda |
| `k` | `NULL` | Infer K from the number of Lambda entries |
| `target_node` | `1` | Node evaluated for Theorem 1 |
| `design_seed` | `238687` | Seed used to construct the latent positions |
| `replication_seed_base` | `2000` | Replication r uses seed `2000 + r` |
| `n_cores` | `7` | Maximum parallel workers |
| `output_dir` | `NULL` | Automatically organize results by setting and model |

Within each setting, the latent positions and Lambda are fixed while the
network observations are generated independently across replications.
Changing `replication_seed_base` changes the observations; changing
`design_seed` changes the latent design. With different n, the same design
seed need not give identical coordinates for node 1. 

The current simulation optimization defaults are:

| Model | Base step-size coefficient (`learning_rate`) | `max_iter` | `tolerance` |
| --- | ---: | ---: | ---: |
| Gaussian | 0.1 | 10000 | 0.001 |
| Poisson | 0.02 | 60000 | 0.001 |
| Logit | 1 | 40000 | 0.001 |
| Probit | 1 | 40000 | 0.001 |

`learning_rate` is a base coefficient. The matrix updates for U and Lambda
multiply their gradients by `learning_rate / n` and `learning_rate / n^2`,
respectively, before projection onto the constraint sets. For example, with
`learning_rate = 1` and `n = 1000`, these coefficients are 0.001 and 0.000001.
The stopping criterion uses the equivalent step `learning_rate / (2*n^2)`
for the scaled `lambda_v` coordinates.

Gaussian uses known noise standard deviation `sigma = 1`. The bounds
`M_bound_U` and `M_bound_Lambda` are both 20. The current tolerance is 0.001
for every model; historical paper runs may use different tolerances. To match
a saved experiment, use its original code and `*_config.rds` settings.
`max_iter` counts optimization steps for **one network**; `n_sim` counts
networks. Models run sequentially, using up to seven workers within each model.

## 2. Reproduce the experiments

Run the examples from the repository root in RStudio. Install any missing
packages:

```r
install.packages(c("rmarkdown", "knitr", "doParallel", "foreach"))
```

To run the current YAML setting:

```r
rmarkdown::render("simulations/node_inference.Rmd")
```

To run all six combinations of sample size and Lambda, with all four models:

```r
lambda_settings <- c(equal = "0.5,0.5", mixed = "0.7,-0.5")
for (n_value in c(250L, 500L, 1000L)) {
  for (setting in names(lambda_settings)) {
    rmarkdown::render(
      "simulations/node_inference.Rmd",
      params = list(
        model = "all", run_simulation = TRUE,
        n = n_value, n_sim = 1000L, k = NULL,
        lambda_values = unname(lambda_settings[setting]),
        target_node = 1L, design_seed = 238687L,
        replication_seed_base = 2000L, n_cores = 7L,
        output_dir = NULL
      ),
      output_file = paste0("node_n", n_value, "_", setting, ".html")
    )
  }
}
```

These commands run the configured experiments. Exact reproduction of saved
paper numbers also requires the corresponding code version and numerical
environment; retain the original results when redrawing published figures.

Results are saved under `data/node1/`, in subfolders for each setting and
model. Keep each `*_results.csv` with its matching `*_config.rds`. Each CSV
row records one attempted replication, including its status. The run also
saves the latent truth, execution settings, complete R results, and basic QQ
plots. Use a new `output_dir` to preserve a separate run.

## 3. Reproduce figures and numerical summaries

### Figures from saved results

Once simulation results are available, run:

```r
rmarkdown::render(
  "plots/plot_figures.Rmd",
  params = list(
    data_dir = "../data/node1",
    output_dir = "../data/figures"
  )
)
```

This creates componentwise normal QQ plots, joint chi-square QQ plots, and
sample-size comparisons for the available settings. Figures are saved as PDF
and PNG. It reads existing results and does not rerun the simulations.
Paths in `params` are relative to the Rmd file: here `../data/node1` refers
to `data/node1` at the repository root. Absolute paths also work.

Keep one consistent experiment collection in `data_dir`; do not mix duplicate
runs or different seeds. Plots retain successful replications with finite
values for every displayed statistic. The input files and retained rows are
recorded alongside the figures.

For all-node boxplots and heatmaps, also supply
`all_nodes_file = "../data/all_nodes/your_run"` in `params`. That directory
must contain complete native `component_summary.csv` files and matching
`config.rds` files. These summaries are generated separately; the single-node
notebook does not produce them. Leave this option unset for single-node QQ plots.

### Numbers for summary tables

Tables report each standardized component separately as **mean (sample
variance)**. Use `var()`, not `sd()`; notebook QQ titles report standard
deviations. The following summarizes one result CSV using the same row
selection as the plotting notebook:

```r
csv <- "path/to/your_results.csv"  # Replace with an actual result file.
d <- read.csv(csv)
components <- grep("^(stat_U_|raw_lambda_|debiased_lambda_)",
                   names(d), value = TRUE)
required <- c(components, "chi_sq_U", "chi_sq_debiased_lambda")
keep <- d$convergence_code == 0 &
  rowSums(!is.finite(as.matrix(d[required]))) == 0
used <- d[which(keep), , drop = FALSE]
stopifnot(nrow(used) >= 2)
summary_table <- data.frame(
  component = components,
  replications = nrow(used),
  mean = vapply(used[components], mean, numeric(1)),
  variance = vapply(used[components], var, numeric(1))
)
summary_table
```

Repeat for each model, n, and Lambda setting. The `raw_lambda_*` and
`debiased_lambda_*` columns contain standardized errors, not parameter
estimates. This produces the component summaries; the plotting notebook does
not export the paper's formatted LaTeX tables.

### Choose HTML or PDF

Both notebooks default to HTML. Use **Knit to PDF** in RStudio, or add
`output_format = "pdf_document"` to `render()`; PDF output requires XeLaTeX.
For a PDF of existing figures, use this option with `plots/plot_figures.Rmd`.
To produce a simulation truth-check PDF without running network replications:

```r
rmarkdown::render(
  "simulations/node_inference.Rmd",
  output_format = "pdf_document",
  params = list(run_simulation = FALSE)
)
```

Changing the report format executes the notebook with the supplied parameters.
`run_simulation = FALSE` only constructs and checks the latent truth; it does
not load previous simulation results.

## 4. Explore other parameter values

Override parameters through `params` in `render()`, or edit them in the
notebook's YAML header. Use `NULL` in R and `null` in YAML for an automatic or
model-default value. For example, a small trial is:

```r
rmarkdown::render(
  "simulations/node_inference.Rmd",
  params = list(
    model = "gaussian", n = 50L, n_sim = 5L,
    lambda_values = "0.5,0.5", k = NULL,
    target_node = 1L, n_cores = 1L, run_simulation = TRUE,
    output_dir = "../data/trial/gaussian"
  ),
  output_file = "trial_gaussian.html"
)
```

This trial checks execution; five replications cannot assess distributional
accuracy. Increase `n_sim` for a substantive simulation comparison and inspect
convergence and truth diagnostics before interpreting the results.

| To change… | Set… |
| --- | --- |
| Observation model | `model = "gaussian"`, `"poisson"`, `"logit"`, or `"probit"` |
| Network size or repetition count | `n` or `n_sim`, respectively |
| Lambda | For example, `lambda_values = "0.7,-0.5"` |
| Latent dimension | For example, `k = NULL, lambda_values = "0.5,0.4,-0.3"` gives K = 3 |
| Node for Theorem 1 | `target_node` to an integer between 1 and n |
| Latent design or network draws | `design_seed` or `replication_seed_base`, respectively |
| Gaussian noise level | `sigma` to a positive standard deviation |
| Available computing resources | `n_cores` from 1 to 7 |

If K is explicit, it must match the number of Lambda entries, with K < n.
All K node components and K(K+1)/2 Lambda components adapt automatically.
Keep the optimizer defaults initially; reducing `max_iter` can produce
nonconverged fits rather than a valid smaller experiment. Save each new
experiment in a separate directory, then point the plotting notebook's
`data_dir` to that directory.
