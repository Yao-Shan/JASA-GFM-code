# Simulations and saved-result figures

## Small trial

Use `simulations/node_inference.Rmd` for one node, or `all_nodes.Rmd` for all nodes.
Both default to `run_simulation: false`, n=60, 3 replications and a single in-process
worker; rendering first checks truth and configuration. Enable `run_simulation`
explicitly to run a trial. Global debiasing is enabled in the single-node notebook
and disabled in the all-node notebook. `n_sim` counts repeated networks, while
`max_iter` limits optimization updates within one network.

```r
rmarkdown::render("simulations/node_inference.Rmd")
# Explicit trial only:
rmarkdown::render("simulations/node_inference.Rmd",
  params = list(run_simulation = TRUE, model = "gaussian", n_sim = 2))
```

New runs reserve unique directories under `data/baseline/alpha_l2/`. Existing run
identifiers cannot be overwritten. Failures and nonconvergence remain in the saved
records; networks are not replaced or selected according to their z statistics.

## Completed-study settings

`simulations/paper_settings.R` returns explicit settings for a selected G1-G5 group
and one Lambda scenario. Sourcing it starts no run. All groups use node1, K2,
alpha817285, design238687, four models, L2-normalized truth alpha, and 1000 attempts
per model/scenario. Model learning rates and iteration limits are Gaussian .1/10000,
Logit and Probit 1/40000, and Poisson .05/80000; tolerance .05, backtracking FALSE.
The full-matrix alpha initializer, bounds20 and probability epsilon .001 are retained.
See `paper_settings.csv` and `docs/seeds.md` for the group mapping.
`docs/study_settings_audit.csv` verifies all14 configuration fields for all40
group/scenario/model combinations against the original saved study settings,
including bounds20, probability epsilon .001 and the n-specific resource allowance.
G2 originally contains other intercept seeds; this package selects817285 explicitly.

```r
source("R/gfm_baseline.R")
source("simulations/paper_settings.R")
args <- gfb_paper_settings("G1", "equal")
str(args)  # inspections only; this does not fit networks
# To execute this full scenario deliberately:
# args$output_dir <- "data/paper/G1_equal"
# run <- do.call(gfb_run_simulation, args)
```

The portable notebook driver computes T1 and T2 from the same fit in one task.
The historical production study scheduled stages separately and preserved fits
between them; the pending-stage controllers, STOP files and machine-specific owner
locks are not shipped as new-run entry points. No existing study is restarted by
this package. For delayed T2 use the saved fit as described in `user_guide.md`.
The notebook driver's preserved limits are at most7 fitting workers and2 inference
slots. They are not a reproduction of historical five-worker T2 scheduling.

At n=2000 each dense T2 calculation can consume several GiB. Study settings retain
6144 MiB per-task resource allowance; the portable helper defaults to one worker.
Worker count can be adjusted explicitly without changing truth or statistical
controls. The published notebooks use ordinary R and do not require Rcpp, a C++
toolchain or a local runtime folder. For multiple workers, select
`worker_backend: "psock"` and set `n_cores` explicitly; the in-process backend
uses one worker.

## Plot saved results

`plots/plot_figures.Rmd` exports separate QQ PDFs and overlaid Lambda histograms.
Its defaults disable plotting and select QQ plots across all saved sample sizes.
Theorem 2 QQ plots show debiased statistics only, without reference bands.
Latin Modern typography requires `pdflatex` and `lmodern`; use `font = "serif"`
for base-R text without TeX.

Run from `with_intercept/`; paths inside `params` are relative to the notebook:

```r
rmarkdown::render("plots/plot_figures.Rmd", params = list(
  run_plots = TRUE,
  results_file = "../data/baseline/alpha_l2/MY_RUN/results.rds",
  n = 1000,
  kinds = c("qq", "hist_estimates", "hist_errors")))
# Historical production exports instead:
# params = list(run_plots = TRUE, data_dir = "../data/MY_EXPORT",
#               models = c("logit", "probit"), n = 1000)
```

Supply either a native `gfm-baseline-v1` results file or one historical experiment
directory containing `theorem1_records.csv` and `theorem2_records.csv`.
Select plots with `kinds`: `qq`, `hist_estimates`, or `hist_errors`.
For node-only runs without global debiasing, select `kinds = "qq"`.
Change `n`, `models`, or `nodes` to choose saved cases; `n: null` keeps all sizes.
`histogram_bins` controls the common bins. Each export creates a new subfolder
under `output_dir`, recording input hashes and retained rows without rerunning
simulations or replacing earlier figures.

The two histogram types compare raw/debiased **aligned matrix estimates** and
those same estimates minus their true entries. The saved coordinate frame is
`Lambda_q`; matrix entry indices are a,b. Diagonal entries of `Lambda_v` are
multiplied by sqrt(2), while off-diagonal entries are unchanged. Standardized
QQ statistics are not converted into parameter estimates. Each histogram pair
uses the same successful, finite replications and bin grid. Density heights are
identical between the estimate and error views; the latter shifts the x-axis by
the true entry.

`plots/plot_individual.R` defines `gfb_read_plot_data()` and `gfb_plot_saved()` for
use without knitting. Sourcing it loads no model or simulation code.

Saved outputs and reports are withheld pending review. Supply saved
per-replication records or run the simulation code to draw QQ plots and
histograms; means and variances alone cannot reconstruct them. To plot a
separately supplied summary CSV:

```sh
Rscript plots/plot_summary.R data/summary.csv my_summary.pdf
```

This draws means and sample variances directly from the supplied summary.
