GFM estimation for one observed network
======================================

This directory provides functions for fitting the symmetric generalized factor
model to an existing dataset and obtaining a bias-corrected Lambda estimate.
The original simulation Rmd files remain unchanged. This interface does not
generate data, repeat a Monte Carlo experiment, choose a seed, or write files.
It uses base R; no package installation is required.

Quick start
-----------
Run these commands from the project root, which contains practitioner/:

    source("practitioner/gfm_estimation.R")
    A <- gfm_read_matrix("network.csv")
    fit <- gfm_fit(A, K = 2, model = "logit", debias = TRUE)

    fit$diagnostics
    fit$debiasing
    fit$Lambda_hat
    fit$Lambda_debiased

Here network.csv is your observed n-by-n network matrix, with no header and no
row labels. It is NOT a CSV of simulation replications or summary statistics.
Choose K for your analysis; the function does not select the latent dimension.

For a CSV whose first row contains column names and first column contains node
names, use:

    A <- gfm_read_matrix("network_with_labels.csv", node_names = TRUE)

The row and column labels must describe the same nodes in the same order.
The reader and fitter do not automatically reorder nodes or fill missing edges.

To fit another observation model:

    fit_probit <- gfm_fit(A, K = 2, model = "probit")
    fit_counts <- gfm_fit(count_matrix, K = 3, model = "poisson")
    fit_normal <- gfm_fit(weight_matrix, K = 2, model = "gaussian", sigma = 1)

Use the matrix appropriate to each model. Gaussian sigma is the known noise
standard deviation, must be supplied explicitly, and is not estimated. Its
default is NULL so it cannot accidentally be treated as known and equal to 1.

Data and model
--------------
For i < j the model has link parameter theta_ij = u_i' Lambda u_j, with
U' U = n I_K and symmetric Lambda. The conditional response means are:

    logit:     1 / (1 + exp(-theta_ij))
    probit:    Phi(theta_ij)
    gaussian:  theta_ij, with conditional variance sigma^2
    poisson:   exp(theta_ij)

A must be a numeric, square, symmetric matrix with complete finite observations
off the diagonal. Both triangles must be supplied. Logit and probit require
binary 0/1 observations; Poisson requires nonnegative integer counts; Gaussian
accepts real observations. The diagonal is ignored and internally set to zero.
Its entries may be zero, nonzero, or NA; they never enter initialization or fitting.
If row and column labels are both present, they must agree in the same order.
The number of nodes n is inferred from A, and 1 <= K < n is required.
Symmetry is checked exactly; the reader never averages mismatched edges.

This implementation assumes a complete undirected network and conditionally
independent edges. It does not implement missing-edge handling, directed data,
self-loop inference, covariates, or a separate intercept/degree-effect term.
Check that these model assumptions suit the dataset before interpreting a fit.

Functions
---------
gfm_fit(A, K, model, sigma = NULL, control = list(), debias = TRUE)
    Fits one network using spectral initialization and projected gradient
    descent. Model names include logit, probit, gaussian, and poisson.
    K is adaptive: Lambda is K-by-K and lambda_v has K(K+1)/2 entries.

gfm_predict(fit, type = "response", debiased = FALSE)
    Returns a symmetric n-by-n matrix in the original node order. The response
    type returns conditional means; type = "link" returns theta_ij. The
    diagonal is NA because self-loops are not modeled. debiased = TRUE uses
    Lambda_debiased and requires a successful correction; it leaves U_hat
    unchanged. These are fitted values for the observed nodes, not predictions
    for new nodes or uncertainty intervals. Both prediction types require a
    converged fit. Correcting Lambda does not establish unbiased edge means.

gfm_read_matrix(file, node_names = FALSE)
    Reads a headerless numeric matrix CSV. With node_names = TRUE, the first
    row and first column instead supply column and row labels.

Optimizer controls
------------------
The defaults match the original simulation optimizer:

    Model       max_iter   learning_rate   tolerance
    logit          40000            1          0.1
    probit         40000            1          0.001
    gaussian       10000            0.1        0.1
    poisson        60000            0.02       0.0001

Additional defaults are M_bound_U = 20, M_bound_Lambda = 20,
probability_eps = 0.001, and max_debias_memory_mb = 2048.
Every control must be one positive finite number. max_iter must be an integer
no larger than R's .Machine$integer.max. probability_eps must be less than 0.5.
Control names must be unique and recognized; misspelled names raise an error.

    fit <- gfm_fit(A, K = 2, model = "logit", control = list(
      max_iter = 60000, learning_rate = 0.5, tolerance = 0.01
    ))

max_iter limits optimization updates on this ONE dataset; it does not specify
Monte Carlo repetitions. learning_rate controls the fixed update steps.
tolerance is the threshold for the scaled projected-gradient residual, not a
confidence level or a relative change in the objective.

M_bound_U bounds each row norm of U and must be at least sqrt(K).
M_bound_Lambda bounds the absolute eigenvalues of Lambda. The joint U
projection is certified only when its polar solution satisfies the row bound;
an active bound that violates this certificate causes an explicit error.
The implementation does not silently replace this projection with row clipping.
probability_eps controls the Bernoulli spectral-initialization clipping, not a
modification of the observed data or fitted likelihood.

Results and status
------------------
The returned list includes:

    U_hat, Lambda_hat, lambda_v_hat
        Raw fitted estimates in the fitted latent coordinate frame. These are
        parameter estimates, not standardized errors from simulation CSVs.
    Lambda_debiased, lambda_v_debiased
        Bias-corrected Lambda estimates; NULL if correction is unavailable.
    correction
        correction_likelihood, correction_penalty, and correction_total;
        all three are vectors in the lambda_v convention below.
    efficient_information
        Information used by the correction, when available.
    diagnostics
        Optimizer convergence and numerical diagnostic information.
    debiasing$status, debiasing$message
        Whether correction succeeded and, if necessary, why it did not.
    config, node_labels
        Resolved model/control settings and retained node identifiers.

Always inspect diagnostics and debiasing before using corrected estimates:

    if (identical(fit$debiasing$status, "success")) {
      corrected_means <- gfm_predict(fit, debiased = TRUE)
    }

Input errors and numerical fitting failures stop with an error. If the
optimizer does not converge, the raw fit is returned with a warning and
debiasing$status is "skipped". If correction fails, including an information,
dimension, or memory check, the raw fit is retained with a warning and status
"failed". Successful correction gives "success". Setting debias = FALSE gives
"not_requested", including when the raw optimizer has not converged.
diagnostics$convergence_code is 0 for convergence and 1 for the iteration
limit; a non-finite objective or residual is a fitting error rather than a
usable returned fit. Do not treat a raw fit returned after a warning as a
successfully corrected estimate; read the accompanying message.

If a fit reaches max_iter, inspect scaled_gradient and the two residuals in
diagnostics. A smaller learning_rate can help when fixed-step updates oscillate;
a larger max_iter can help when progress is slow. Refit the same observed A
with explicitly chosen controls. Do not loosen tolerance merely to force a
successful status. Convergence is a stationarity check, not proof that this
nonconvex problem has reached a global optimum. An information/rank failure
can indicate that the chosen model or K is not identifiable from this network.

Coordinate convention and interpretation
---------------------------------------
The model is unchanged by an orthogonal rotation of U accompanied by the
corresponding conjugation of Lambda. With real data, the true U is unavailable
for the alignment used in the simulation study. This interface uses the fitted
U_hat as the plug-in constraint anchor and reports correction in that frame.
It does not reproduce the simulation's unobservable truth-aligned errors and
does not provide coordinatewise confidence intervals with a claimed theorem
guarantee for this plug-in frame. U_hat itself is not bias-corrected.

lambda_v stores the K diagonal elements of Lambda divided by sqrt(2), followed
by the strict upper-triangular elements in row order: (1,2), (1,3), ..., then
(2,3), ..., (K-1,K). For K = 2 it is
(Lambda_11/sqrt(2), Lambda_22/sqrt(2), Lambda_12). The matrix-valued outputs
convert this convention back to ordinary Lambda entries.
The reported total correction is an additive adjustment:
lambda_v_debiased = lambda_v_hat + total correction.
The corrected Lambda is not reprojected onto the optimizer's spectral bound;
its eigenvalues can exceed M_bound_Lambda after the adjustment.

Computation, saving, and checks
------------------------------
The implementation uses dense matrices. Fitting includes full spectral
decompositions, and correction forms matrices of dimension nK-by-nK with
cubic-cost linear algebra. The correction guard estimates workspace as sixteen
dense (nK)-by-(nK) double arrays and compares it with max_debias_memory_mb.
Despite the control's name, this limit is measured in MiB (1024^2 bytes).
This estimate is not a guaranteed peak-memory bound. Check available resources
before raising the limit; this is not a distributed or large-network solver.
Use debias = FALSE when only the raw fit is needed.

Sourcing the file only defines functions. Fitting does not change the random
seed, launch a simulation, create files, or start parallel Monte Carlo workers.
Saving is an explicit, optional action:

    saveRDS(fit, "network_fit.rds")

Run the supplied deterministic checks from the project root with:

    Rscript --vanilla practitioner/tests/test_gfm.R
