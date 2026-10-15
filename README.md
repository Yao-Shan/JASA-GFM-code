# JASA GFM code

R implementations of a symmetric generalized factor model for Gaussian,
Poisson, Bernoulli-logit, and Bernoulli-probit network observations.

The repository provides an interface for fitting **one existing network**, a
separate Monte Carlo simulation notebook, and a plotting notebook for saved
simulation results.

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
│   ├── user_guide.md          # Inputs, controls, results, and interpretation
│   └── simulation_guide.md    # Paper settings, small trial runs, and plotting
├── plots/
│   └── plot_figures.Rmd       # Figures from existing results; no fitting
├── simulations/
│   └── node_inference.Rmd     # Data generation and Monte Carlo analysis
├── tests/
│   ├── test_gfm.R             # Estimation and formula regression checks
│   └── test_simulation_config.R # Defaults, dimensions, and output paths
└── README.md
```

For an existing dataset, start with [`R/gfm_estimation.R`](R/gfm_estimation.R)
and the [user guide](docs/user_guide.md). To run the simulation study, use
[`simulations/node_inference.Rmd`](simulations/node_inference.Rmd) and the
[simulation guide](docs/simulation_guide.md). To redraw saved results, use
[`plots/plot_figures.Rmd`](plots/plot_figures.Rmd).

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
Rscript --vanilla tests/test_simulation_config.R
```

The estimation suite covers all four likelihoods, agreement with
the simulation formulas under the same coordinate anchor, adaptive dimensions,
rotation invariance, input validation, failure handling, and predictions.
Fixtures are used only inside the test script. The suite extracts selected
function definitions from the notebook without running its simulation pipeline.
The configuration suite also checks the YAML parameters and requires `knitr`;
it does not start simulations. The estimation suite requires only base R.

## Simulate networks or redraw saved results

Open `simulations/node_inference.Rmd` in RStudio to reproduce or adapt the simulation
workflow. It uses `rmarkdown`, `knitr`, `doParallel`, and `foreach`, which are
not required by the practitioner interface or the estimation test suite.

The defaults retain the paper's main setting: all four models, `n = 1000`,
1,000 network replications, node 1, Lambda diagonal `(0.5, 0.5)`, and up to seven
workers. K is inferred from the supplied Lambda entries. Knitting with
`run_simulation: true` starts that experiment; set it to `false` to construct
and check the truth without running network replications.

The [simulation guide](docs/simulation_guide.md) gives the two Lambda settings
at `n = 250, 500, 1000`, a small trial configuration, seed provenance, and the
difference between network repetitions and optimizer iterations. It also
explains automatic output paths and how to retain separate runs.

To redraw existing results without fitting models:

```r
rmarkdown::render("plots/plot_figures.Rmd")
```

Both notebooks default to HTML and also support **Knit to PDF** or
`output_format = "pdf_document"`. PDF rendering requires XeLaTeX and retains
the generated `.tex` source. PDF reports hide code; HTML reports offer
collapsible code. The plotting report embeds vector PDF figures, and separate
PDF/PNG figure exports remain available. Choosing another format reruns the
notebook with its current parameters; see the
[report-format examples](docs/simulation_guide.md#choose-html-or-pdf) for a
simulation truth-check PDF and a PDF made from existing results.

The plotting notebook reads native result CSVs and their matching configuration
RDS files from `data/node1` by default. It writes open-circle QQ plots and
available sample-size comparisons to `data/figures`, with records of the input
files and retained replications. Each component is displayed separately, and
dimensions adapt to K. Optional all-node figures require separate all-node
summaries and metadata; the nodewise simulation notebook does not generate
those summaries.

Simulation results, all-node summaries, and analysis datasets are generated or
supplied separately; they are not bundled in this code repository. The code
alone is therefore not a bundle of every saved paper figure and dataset.
