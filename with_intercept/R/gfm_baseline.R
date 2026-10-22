# Baseline symmetric GFM: alpha_i + alpha_j + u_i' Lambda u_j.
# Source this entry point to load functions without running a simulation.
.gfb_source_candidates <- vapply(sys.frames(), function(frame) {
  file <- frame$ofile
  if (is.character(file) && length(file) == 1L && !is.na(file) &&
      basename(file) == "gfm_baseline.R" && file.exists(file)) file else NA_character_
}, character(1))
.gfb_source_candidates <- .gfb_source_candidates[!is.na(.gfb_source_candidates)]
if (!length(.gfb_source_candidates)) stop("Load this file with source().")
.gfb_entry_file <- tail(.gfb_source_candidates, 1L)
.gfb_root <- dirname(normalizePath(.gfb_entry_file, mustWork = TRUE))
source(file.path(.gfb_root, "gfm_baseline_utilities.R"), local = TRUE)
rm(.gfb_entry_file, .gfb_source_candidates)

.gfb_scalar <- function(x, name, lower = -Inf, integer = FALSE, strict = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      (if (strict) x <= lower else x < lower) ||
      (integer && x != floor(x))) stop("Invalid ", name, ".")
  invisible(x)
}

# Validate model and optimizer settings.
gfb_config <- function(n, K = 2L, model = "gaussian", sigma = 1,
                       control = list()) {
  .gfb_scalar(n, "n", 3, TRUE)
  .gfb_scalar(K, "K", 1, TRUE)
  if (K >= n) stop("Centered coordinates require K < n.")
  if (!is.character(model) || length(model) != 1L || is.na(model))
    stop("model must name one observation model.")
  model <- .gfb_normalize_model_name(model)
  .gfb_scalar(sigma, "sigma", 0, strict = TRUE)
  defaults <- list(
    learning_rate = switch(model, gaussian = .1, poisson = .05, 1),
    max_iter = switch(model, gaussian = 10000L, poisson = 80000L, 40000L),
    tolerance = 1e-3, M_bound_U = 20, M_bound_Lambda = 20,
    M_bound_alpha = 20, probability_eps = 1e-3,
    max_debias_memory_mb = 2048, backtracking = FALSE,
    max_backtracks = 30L
  )
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || any(names(control) == "") ||
       anyDuplicated(names(control)) || any(!names(control) %in% names(defaults)))))
    stop("control must be a uniquely named list of supported controls.")
  for (key in names(control)) defaults[[key]] <- control[[key]]
  for (key in c("learning_rate", "tolerance", "M_bound_U", "M_bound_Lambda",
                "M_bound_alpha", "max_debias_memory_mb"))
    .gfb_scalar(defaults[[key]], key, 0, strict = TRUE)
  if (defaults$M_bound_U < sqrt(K)) stop("M_bound_U must be at least sqrt(K).")
  .gfb_scalar(defaults$max_iter, "max_iter", 1, TRUE)
  .gfb_scalar(defaults$max_backtracks, "max_backtracks", 0, TRUE)
  .gfb_scalar(defaults$probability_eps, "probability_eps", 0, strict = TRUE)
  if (defaults$probability_eps >= .5) stop("probability_eps must be below 0.5.")
  if (!is.logical(defaults$backtracking) || length(defaults$backtracking) != 1L ||
      is.na(defaults$backtracking)) stop("backtracking must be TRUE or FALSE.")
  c(list(n = as.integer(n), k = as.integer(K), model = model, sigma = sigma), defaults)
}

# Run seeded draws and restore the caller's RNG state.
.gfb_with_seed <- function(seed, expr) {
  .gfb_scalar(seed, "seed", 0, TRUE)
  if (seed > .Machine$integer.max) stop("seed is too large.")
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  eval(substitute(expr), envir = parent.frame())
}

.gfb_theta <- function(U, alpha, Lambda) {
  U %*% Lambda %*% t(U) + outer(alpha, alpha, "+")
}

.gfb_safe_sqrt_psd <- function(A) {
  e <- eigen(.gfb_symmetrize(A), symmetric = TRUE)
  if (any(!is.finite(e$values)) || min(e$values) < -1e-9 * max(1, abs(e$values)))
    stop("Matrix is not positive semidefinite.")
  .gfb_symmetrize(sweep(e$vectors, 2L, sqrt(pmax(e$values, 0)), "*") %*% t(e$vectors))
}

# Optionally normalize the sampled alpha vector to unit L2 norm.
# Fitted alpha values are not normalized.
.gfb_prepare_alpha <- function(alpha_raw, n, alpha_normalization = c("l2", "none")) {
  alpha_normalization <- match.arg(alpha_normalization)
  .gfb_scalar(n, "alpha dimension", 1, TRUE)
  if (!is.numeric(alpha_raw) || !is.null(dim(alpha_raw)) || length(alpha_raw) != n ||
      any(!is.finite(alpha_raw))) stop("alpha_raw must be a finite numeric vector of length n.")
  scale <- max(abs(alpha_raw))
  # Scale first to compute the norm stably.
  raw_norm <- if (scale == 0) 0 else scale * sqrt(sum((alpha_raw / scale)^2))
  if (!is.finite(raw_norm)) stop("The Euclidean norm of alpha_raw must be finite.")
  if (alpha_normalization == "l2" && raw_norm <= 0)
    stop("L2 normalization requires a strictly positive Euclidean norm of alpha_raw.")
  denominator <- if (alpha_normalization == "l2") raw_norm else 1
  alpha <- if (alpha_normalization == "l2") alpha_raw / denominator else alpha_raw
  list(alpha = alpha, raw_norm = raw_norm,
       metadata = list(method = alpha_normalization, dimension = as.integer(n),
         denominator = denominator, raw_l2_norm = raw_norm,
         l2_norm = if (alpha_normalization == "l2") sqrt(sum(alpha^2)) else raw_norm,
         formula = if (alpha_normalization == "l2") "alpha_raw / sqrt(sum(alpha_raw^2))" else "alpha_raw"))
}

# Generate fixed U_star, alpha_star and Lambda_star for a study.
gfb_generate_truth <- function(n = 60L, K = 2L, lambda_values = c(.5, .5),
                               design_seed = 238687L, alpha_seed = 817263L,
                               alpha_min = -2, alpha_max = 0,
                               alpha_normalization = c("l2", "none")) {
  alpha_normalization <- match.arg(alpha_normalization)
  config <- gfb_config(n, K)
  n <- config$n; K <- config$k
  if (!is.numeric(lambda_values) || length(lambda_values) != K ||
      any(!is.finite(lambda_values)) || any(lambda_values == 0))
    stop("lambda_values must contain K finite nonzero values.")
  .gfb_scalar(alpha_min, "alpha_min")
  .gfb_scalar(alpha_max, "alpha_max")
  if (alpha_max <= alpha_min) stop("alpha_max must exceed alpha_min.")
  .gfb_scalar(design_seed, "design_seed", 0, TRUE)
  .gfb_scalar(alpha_seed, "alpha_seed", 0, TRUE)
  if (design_seed == alpha_seed) stop("Use distinct design_seed and alpha_seed.")
  sampled <- .gfb_with_seed(design_seed, {
    Z <- matrix(0, n, K); filled <- 0L; proposed <- accepted_count <- 0
    batch_size <- as.integer(min(50000, n * 2^min(K, 15)))
    while (filled < n) {
      candidates <- matrix(runif(batch_size * K, -1, 1), batch_size, K)
      accepted <- candidates[rowSums(candidates^2) < 1, , drop = FALSE]
      proposed <- proposed + batch_size
      accepted_count <- accepted_count + nrow(accepted)
      take <- min(n - filled, nrow(accepted))
      if (take > 0L) Z[filled + seq_len(take), ] <- accepted[seq_len(take), , drop = FALSE]
      filled <- filled + take
      if (proposed > 2e7) stop("Ball rejection sampler exceeded its proposal budget.")
    }
    list(Z = Z, proposed = proposed, accepted = accepted_count)
  })
  Zc <- sweep(sampled$Z, 2L, colMeans(sampled$Z), "-")
  e <- eigen(crossprod(Zc) / n, symmetric = TRUE)
  if (min(e$values) <= 1e-10 * max(e$values)) stop("Centered second moment is singular.")
  invsqrt <- sweep(e$vectors, 2L, 1 / sqrt(e$values), "*") %*% t(e$vectors)
  U <- Zc %*% invsqrt
  alpha_raw <- .gfb_with_seed(alpha_seed, runif(n, alpha_min, alpha_max))
  alpha_prepared <- .gfb_prepare_alpha(alpha_raw, n, alpha_normalization)
  alpha <- alpha_prepared$alpha
  raw_distribution <- list(name = "uniform", min = alpha_min, max = alpha_max)
  distribution <- if (alpha_normalization == "none") raw_distribution else
    list(name = "l2_normalized_uniform_vector", raw_distribution = raw_distribution,
         normalization = "l2", dimension = n, target_l2_norm = 1, iid = FALSE)
  Lambda <- diag(lambda_values, K, K)
  Theta <- .gfb_theta(U, alpha, Lambda)
  list(U_star = U, alpha_star = alpha, Lambda_star = Lambda, Theta_star = Theta,
       seeds = list(design_seed = design_seed, alpha_seed = alpha_seed),
       alpha_raw = alpha_raw, alpha_raw_norm = alpha_prepared$raw_norm,
       alpha_raw_distribution = raw_distribution,
       alpha_normalization = alpha_normalization,
       alpha_normalization_metadata = alpha_prepared$metadata,
       alpha_distribution = distribution,
       diagnostics = list(center_error = max(abs(colMeans(U))),
         gram_error = max(abs(crossprod(U) / n - diag(K))),
         max_row_norm = max(sqrt(rowSums(U^2))),
         rejection_acceptance_rate = sampled$accepted / sampled$proposed,
         theta_range = range(Theta[upper.tri(Theta)])))
}

# Sample one network conditional on the fixed truth.
gfb_generate_network <- function(truth, model = "gaussian", sigma = 1, seed = 2001L) {
  if (!is.list(truth) || !is.matrix(truth$U_star)) stop("Invalid truth object.")
  config <- gfb_config(nrow(truth$U_star), ncol(truth$U_star), model, sigma)
  Theta <- .gfb_theta(truth$U_star, truth$alpha_star, truth$Lambda_star)
  if (any(!is.finite(Theta))) stop("Non-finite truth predictor.")
  upper <- upper.tri(Theta)
  th <- Theta[upper]
  values <- .gfb_with_seed(seed, {
    switch(config$model,
      gaussian = rnorm(length(th), th, config$sigma),
      poisson = { intensity <- exp(th)
        if (any(!is.finite(intensity))) stop("Non-finite Poisson means.")
        rpois(length(th), intensity) },
      bernoulli_logit = rbinom(length(th), 1L, plogis(th)),
      bernoulli_probit = rbinom(length(th), 1L, pnorm(th)))
  })
  A <- matrix(0, config$n, config$n)
  A[upper] <- values
  A <- A + t(A)
  A
}

# Project onto centered U with U'U = nI.
# Report a violated row bound; clipping would break U'U = nI.
.gfb_project_U <- function(Z, bound = 20) {
  if (!is.matrix(Z) || !is.numeric(Z) || any(!is.finite(Z))) stop("Invalid U update.")
  n <- nrow(Z); K <- ncol(Z)
  .gfb_scalar(bound, "U row bound", sqrt(K))
  if (K < 1L || n <= K) stop("Centered U projection requires n > K >= 1.")
  Zc <- sweep(Z, 2L, colMeans(Z), "-")
  s <- svd(Zc, nu = K, nv = K)
  if (min(s$d) <= 1e-12 * max(1, s$d)) stop("Centered U update is rank deficient.")
  U <- sqrt(n) * s$u %*% t(s$v)
  center_error <- max(abs(colMeans(U)))
  gram_error <- max(abs(crossprod(U) / n - diag(K)))
  max_row_norm <- max(sqrt(rowSums(U^2)))
  objective <- .5 * sum((U - Z)^2)
  lower_bound <- .5 * (sum(Z^2) + n * K - 2 * sqrt(n) * sum(s$d))
  objective_gap <- abs(objective - lower_bound)
  if (center_error > 1e-9 || gram_error > 1e-9 ||
      objective_gap > 1e-8 * max(1, abs(objective), abs(lower_bound)))
    stop("Centered polar projection failed its numerical certificate.")
  if (max_row_norm > bound + 1e-10 * max(1, bound))
    stop("U row bound is active; a certified joint projection is unavailable.")
  list(U = U, center_error = center_error, gram_error = gram_error,
       max_row_norm = max_row_norm, row_margin = bound - max_row_norm,
       objective_gap = objective_gap, method = "centered_polar_exact")
}

.gfb_objective <- function(U, alpha, Lambda, A, config) {
  Theta <- .gfb_theta(U, alpha, Lambda)
  upper <- upper.tri(A)
  .gfb_negative_loglikelihood_edges(Theta[upper], A[upper], config)
}

.gfb_gradient <- function(U, alpha, Lambda, A, config) {
  D <- .gfb_negative_score_matrix(.gfb_theta(U, alpha, Lambda), A, config)
  diag(D) <- 0
  list(U = D %*% U %*% Lambda, alpha = rowSums(D),
       Lambda = .gfb_symmetrize(.5 * crossprod(U, D %*% U)))
}

# Initialize alpha = (n I + 11')^{-1} B 1 using the full reconstructed B,
# including its diagonal; the likelihood uses off-diagonal edges only.
.gfb_additive_initial <- function(B) {
  n <- nrow(B)
  rowSums(B) / n - sum(B) / (2 * n^2)
}

.gfb_initialize <- function(A, config) {
  n <- config$n; K <- config$k
  s <- svd(A)
  # Retain K+2 directions to include the rank-two alpha contribution.
  r <- min(n, K + 2L)
  keep <- .gfb_spectral_keep_indices(s$d, r, if (config$model == "gaussian") 0L else 2L)
  mean_est <- .gfb_symmetrize(s$u[, keep, drop = FALSE] %*%
    (s$d[keep] * t(s$v[, keep, drop = FALSE])))
  B <- switch(config$model,
    gaussian = mean_est,
    poisson = log(pmax(mean_est, .05)),
    bernoulli_logit = qlogis(pmin(pmax(mean_est, config$probability_eps), 1-config$probability_eps)),
    bernoulli_probit = qnorm(pmin(pmax(mean_est, config$probability_eps), 1-config$probability_eps)))
  alpha <- .gfb_additive_initial(B)
  if (any(!is.finite(alpha)) || max(abs(alpha)) > config$M_bound_alpha)
    stop("Full-matrix alpha initializer is non-finite or exceeds M_bound_alpha; it was not clipped.")
  # Remove centering roundoff, then initialize the signed rank-K signal.
  # Keep alpha from the full-matrix normal equations.
  residual <- B - outer(alpha, alpha, "+")
  centered <- sweep(sweep(residual, 1L, rowMeans(residual), "-"),
                    2L, colMeans(residual), "-") + mean(residual)
  e <- eigen(.gfb_symmetrize(centered), symmetric = TRUE)
  idx <- order(abs(e$values), decreasing = TRUE)[seq_len(K)]
  U <- .gfb_project_U(sqrt(n) * e$vectors[, idx, drop = FALSE], config$M_bound_U)$U
  Lambda <- .gfb_proximal_projection_Lambda(diag(e$values[idx] / n, K, K), config$M_bound_Lambda)
  list(U = U, alpha = alpha, Lambda = Lambda)
}

.gfb_step <- function(U, alpha, Lambda, gradient, config, learning_rate) {
  n <- config$n
  projection <- .gfb_project_U(U - learning_rate * gradient$U / n, config$M_bound_U)
  next_alpha <- pmin(pmax(alpha - learning_rate * gradient$alpha / n,
                          -config$M_bound_alpha), config$M_bound_alpha)
  next_Lambda <- .gfb_proximal_projection_Lambda(
    Lambda - learning_rate * gradient$Lambda / n^2, config$M_bound_Lambda)
  # For lambda_v = (diag(Lambda)/sqrt(2), offdiag(Lambda)), use lr/(2 n^2).
  residuals <- c(U = max(abs(projection$U - U)) / (learning_rate / n),
    alpha = max(abs(next_alpha - alpha)) / (learning_rate / n),
    Lambda_v = max(abs(.gfb_get_lambda_v(next_Lambda) - .gfb_get_lambda_v(Lambda))) /
      (learning_rate / (2 * n^2)))
  list(U = projection$U, alpha = next_alpha, Lambda = next_Lambda,
       residuals = residuals, projection = projection)
}

# Fit one observed network and return estimates plus convergence diagnostics.
gfb_fit <- function(A, K = 2L, model = "gaussian", sigma = 1, control = list()) {
  model <- .gfb_normalize_model_name(model)
  A <- .gfb_validate_matrix(A, model)
  config <- gfb_config(nrow(A), K, model, sigma, control)
  current <- .gfb_initialize(A, config)
  U <- current$U; alpha <- current$alpha; Lambda <- current$Lambda
  objective <- .gfb_objective(U, alpha, Lambda, A, config)
  if (!is.finite(objective)) stop("Non-finite initialized likelihood.")
  code <- 1L; n_iter <- 0L; backtracks <- 0L; failure_message <- ""
  accepted_lr <- config$learning_rate
  # Evaluate the stopping residual at the returned parameters.
  for (iteration in seq_len(config$max_iter)) {
    gradient <- .gfb_gradient(U, alpha, Lambda, A, config)
    lr <- config$learning_rate
    accepted <- FALSE
    for (bt in 0:config$max_backtracks) {
      proposal <- tryCatch(.gfb_step(U, alpha, Lambda, gradient, config, lr),
                           error = identity)
      next_objective <- if (inherits(proposal, "error")) Inf else
        .gfb_objective(proposal$U, proposal$alpha, proposal$Lambda, A, config)
      if (inherits(proposal, "error")) failure_message <- conditionMessage(proposal)
      if (is.finite(next_objective) && (!config$backtracking ||
          next_objective <= objective + 1e-12 * max(1, abs(objective)))) {
        accepted <- TRUE
        break
      }
      if (!config$backtracking) break
      lr <- lr / 2
      backtracks <- backtracks + 1L
    }
    if (!accepted) {
      if (!nzchar(failure_message)) failure_message <- if (config$backtracking)
        "No finite decreasing projected step found." else
        "No finite projected step found at the fixed learning rate."
      code <- 2L; break
    }
    accepted_lr <- lr
    if (max(proposal$residuals) <= config$tolerance) { code <- 0L; break }
    U <- proposal$U; alpha <- proposal$alpha; Lambda <- proposal$Lambda
    objective <- next_objective
    n_iter <- iteration
  }
  final_gradient <- .gfb_gradient(U, alpha, Lambda, A, config)
  final_step <- tryCatch(.gfb_step(U, alpha, Lambda, final_gradient, config, accepted_lr),
                         error = identity)
  if (inherits(final_step, "error")) {
    code <- 2L
    failure_message <- conditionMessage(final_step)
    residuals <- c(U = Inf, alpha = Inf, Lambda_v = Inf)
  } else residuals <- final_step$residuals
  if (code != 2L) code <- if (max(residuals) <= config$tolerance) 0L else 1L
  projection <- .gfb_project_U(U, config$M_bound_U)
  one_sided_nodes <- if (model == "gaussian") integer() else {
    low <- which(rowSums(A) == 0)
    high <- if (model %in% c("bernoulli_logit", "bernoulli_probit"))
      which(rowSums(A) == config$n - 1L) else integer()
    sort(unique(c(low, high)))
  }
  diagnostics <- list(convergence_code = code, converged = code == 0L,
    failure_message = if (code == 2L) failure_message else "",
    one_sided_nodes = one_sided_nodes,
    n_iter = n_iter, scaled_gradient = max(residuals),
    projected_residual_U = unname(residuals["U"]),
    projected_residual_alpha = unname(residuals["alpha"]),
    projected_residual_Lambda_v = unname(residuals["Lambda_v"]),
    final_objective = .gfb_objective(U, alpha, Lambda, A, config),
    effective_learning_rate = accepted_lr, backtracking_count = backtracks,
    center_error = projection$center_error, gram_error = projection$gram_error,
    max_row_norm = projection$max_row_norm,
    alpha_bound_margin = config$M_bound_alpha - max(abs(alpha)),
    lambda_bound_margin = config$M_bound_Lambda - max(abs(eigen(Lambda, symmetric=TRUE, only.values=TRUE)$values)))
  structure(list(U_hat = U, alpha_hat = alpha, Lambda_hat = Lambda,
    A = A, config = config, diagnostics = diagnostics), class = "gfb_fit")
}

# Return fitted link values or response means.
gfb_predict <- function(fit, type = c("response", "link")) {
  type <- match.arg(type)
  if (!inherits(fit, "gfb_fit")) stop("fit must be a gfb_fit object.")
  if (fit$diagnostics$convergence_code != 0L) stop("Prediction requires a converged fit.")
  Theta <- .gfb_theta(fit$U_hat, fit$alpha_hat, fit$Lambda_hat)
  out <- if (type == "link") Theta else .gfb_link_probability(Theta, fit$config)
  diag(out) <- NA_real_
  out
}

print.gfb_fit <- function(x, ...) {
  cat("Baseline symmetric GFM:", x$config$model, "n =", x$config$n, "K =", x$config$k, "\n")
  cat("Converged:", x$diagnostics$converged, "; iterations:", x$diagnostics$n_iter,
      "; residual:", format(x$diagnostics$scaled_gradient), "\n")
  invisible(x)
}

# Load base-R helpers into this environment.
for (.gfb_module in c("gfm_baseline_inference.R", "gfm_baseline_simulation.R")) {
  .gfb_module_path <- file.path(.gfb_root, .gfb_module)
  source(.gfb_module_path, local = TRUE)
}
rm(.gfb_module, .gfb_module_path)

# Record the method used by workers and saved results.
.gfb_debias_method <- "pdf_bp_bm"
