# Deterministic regression checks; run with Rscript --vanilla from any directory.
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
# R may encode spaces as ~+~ in its --file command-line argument.
if (length(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(script[[1]]))) else getwd()
source_file <- file.path(root, "R", "gfm_estimation.R")
stopifnot(file.exists(source_file))
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAIL: ", label, call. = FALSE)
  checks <<- checks + 1L
}
equal <- function(x, y, label, tolerance = 1e-10) {
  result <- all.equal(x, y, tolerance = tolerance, check.attributes = FALSE)
  if (!isTRUE(result)) stop("FAIL: ", label, ": ", paste(result, collapse = "; "), call. = FALSE)
  checks <<- checks + 1L
}
error <- function(expr, pattern) {
  caught <- tryCatch({ force(expr); NULL }, error = identity)
  check(inherits(caught, "error") && grepl(pattern, conditionMessage(caught)),
        paste("expected error matching", pattern))
}
warnings_from <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning")
  })
  list(value = value, messages = messages)
}
snapshot <- function() list(kind = RNGkind(), options = options(), search = search(),
  seed = if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) .Random.seed else NULL)
api <- new.env(parent = globalenv())
before <- snapshot()
sys.source(source_file, envir = api)
check(identical(before, snapshot()), "source has no RNG/options/search side effects")
set.seed(712) # Test fixture only; the production source never sets a seed.
before <- snapshot()
sys.source(source_file, envir = api)
check(identical(before, snapshot()), "source preserves an existing random seed")
check(all(vapply(as.list(api), is.function, logical(1))), "source defines functions only")
fit_network <- api$gfm_fit
predict_network <- api$gfm_predict

# Load selected function definitions from R chunks without running the pipeline.
lines <- readLines(file.path(root, "simulations", "node_inference.Rmd"), warn = FALSE)
in_chunk <- FALSE
r_code <- character()
for (line in lines) {
  if (grepl("^```\\{r", line)) { in_chunk <- TRUE; next }
  if (in_chunk && grepl("^```\\s*$", line)) { in_chunk <- FALSE; next }
  if (in_chunk) r_code <- c(r_code, line)
}
reference <- new.env(parent = globalenv())
original_names <- sub("^\\.gfm_", "", ls(api, all.names = TRUE))
for (expr in parse(text = r_code)) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.symbol(expr[[2]]) && is.call(expr[[3]]) &&
      identical(expr[[3]][[1]], as.name("function")) &&
      as.character(expr[[2]]) %in% original_names) eval(expr, reference)
}
check(exists("run_pgd_network", reference, inherits = FALSE), "original PGD loaded safely")

# Check explicit n_drop inverses in the source and notebook helpers.
for (implementation in c("source", "simulation")) {
  pinv <- if (implementation == "source") api$.gfm_safe_pinv_psd else reference$safe_pinv_psd
  check(identical(names(formals(pinv)), c("A", "n_drop")) &&
        identical(formals(pinv)$n_drop, 0L),
        paste(implementation, "pseudoinverse exposes only n_drop"))
  tiny <- pinv(diag(c(2, 1e-16)))
  equal(tiny[2, 2], 1e16,
        paste(implementation, "retains a tiny strictly positive eigenvalue"))
  equal(diag(c(2, 1e-16)) %*% tiny, diag(2),
        paste(implementation, "full-rank inverse identity"))
  equal(pinv(diag(c(4, 0, 2)), n_drop = 1L), diag(c(0.25, 0, 0.5)),
        paste(implementation, "removes one explicitly specified null direction"))
  equal(pinv(matrix(1, 2, 2), n_drop = 1L), matrix(0.25, 2, 2),
        paste(implementation, "rank-one inverse in rotated coordinates"))
  equal(pinv(diag(c(4, -0.01)), n_drop = 1L), diag(c(0.25, 0)),
        paste(implementation, "drops the eigenvalue closest to zero"))
  equal(pinv(diag(c(-2, 0)), n_drop = 2L), matrix(0, 2, 2),
        paste(implementation, "dropping every direction returns zero"))
  error(pinv(diag(c(1, 0))), "Remaining eigenvalues must be strictly positive")
  error(pinv(diag(c(1, 0, 0)), n_drop = 1L),
        "Remaining eigenvalues must be strictly positive")
  error(pinv(diag(c(1, -0.1))), "Remaining eigenvalues must be strictly positive")
  error(pinv(diag(c(0.01, -2)), n_drop = 1L),
        "Remaining eigenvalues must be strictly positive")
  for (invalid_drop in list(-1L, 3L, 0.5, numeric(), c(0L, 1L), Inf, NA_real_)) {
    error(pinv(diag(2), n_drop = invalid_drop), "n_drop must be an integer")
  }
}

# One observed-network fixture, not a Monte Carlo experiment.
n <- 35L; K <- 2L
U <- sqrt(n) * qr.Q(qr(matrix(rnorm(n * K), n, K)))
theta <- U %*% diag(c(0.5, -0.4)) %*% t(U)
A <- matrix(0, n, n)
A[upper.tri(A)] <- rnorm(n * (n - 1) / 2, theta[upper.tri(theta)], 0.5)
A <- A + t(A)
dimnames(A) <- rep(list(paste0("node_", seq_len(n))), 2)
before <- snapshot()
fit <- fit_network(A, K, "gaussian", sigma = 0.5, U_star = U)
check(identical(before, snapshot()), "fit has no RNG/options/search side effects")
check(fit$diagnostics$converged && fit$debiasing$status == "success", "Gaussian fit and correction succeed")
equal(crossprod(fit$U_hat) / n, diag(K), "fitted U normalization")
equal(fit$lambda_v_debiased, fit$lambda_v_hat + fit$correction$correction_total,
      "additive correction sign")
equal(api$.gfm_get_lambda_v(fit$Lambda_debiased), fit$lambda_v_debiased,
      "corrected matrix/vector correspondence")

# Compare initialization, updates and correction across all four likelihoods.
observations <- list(gaussian = A, logit = 1 * (A > 0),
                     probit = 1 * (A > 0), poisson = round(abs(A)))
for (model in names(observations)) {
  B <- observations[[model]]; diag(B) <- 0
  config <- api$.gfm_config(n, K, model, if (model == "gaussian") 0.5 else NULL,
                            list(max_iter = 5L, tolerance = 1e-12))
  equal(api$.gfm_initialize_spectral(B, config), reference$initialize_spectral(B, config),
        paste(model, "spectral initialization parity"))
  equal(api$.gfm_run_pgd_network(B, config), reference$run_pgd_network(B, config),
        paste(model, "five-step PGD parity"))
  modern <- api$.gfm_calculate_debiasing_components(fit$U_hat, fit$Lambda_hat,
                                                  U, config, B)
  old <- reference$calculate_debiasing_components(fit$U_hat, fit$Lambda_hat,
                                                 U, config, B)
  equal(modern, old, paste(model, "correction parity with true anchor"))
}
raw <- reference$run_pgd_network(A, fit$config)
aligned <- reference$rotate_estimator(raw$U_hat, raw$Lambda_hat, U)
equal(aligned$U_q, fit$U_hat, "truth-aligned U matches simulation")
equal(aligned$Lambda_q, fit$Lambda_hat, "truth-aligned Lambda matches simulation")
equal(aligned$Q_transpose, fit$alignment$rotation, "alignment rotation matches simulation")
expected_correction <- reference$calculate_debiasing_components(
  aligned$U_q, aligned$Lambda_q, U, fit$config, A)
equal(fit$correction$correction_total, expected_correction$correction_total,
      "public correction uses true U, matching simulation")
raw_only <- fit_network(A, K, "gaussian", sigma = 0.5, debias = FALSE)
equal(raw_only$U_hat, raw$U_hat, "truth-free raw U matches optimizer")
equal(raw_only$Lambda_hat, raw$Lambda_hat, "truth-free raw Lambda matches optimizer")
equal(predict_network(raw_only), predict_network(fit), "alignment preserves predictions")
check(raw_only$alignment$status == "not_requested", "raw-only fit needs no truth")

# Verify adaptive parameter ordering beyond the K=2 simulation case.
for (k in c(1L, 3L)) {
  L <- matrix(seq_len(k * k), k, k); L <- (L + t(L)) / 2
  v <- api$.gfm_get_lambda_v(L)
  equal(api$.gfm_lambda_matrix(v, k), L, paste("K", k, "parameter round trip"))
  check(length(v) == k * (k + 1) / 2, paste("K", k, "parameter count"))
  x <- seq_len(k) / 3; y <- rev(seq_len(k)) / 4
  equal(sum(v * api$.gfm_lambda_direction_derivative(x, y)),
        as.numeric(t(x) %*% L %*% y), paste("K", k, "predictor parameterization"))
  set.seed(900 + k)
  Uk <- sqrt(24) * qr.Q(qr(matrix(rnorm(24 * k), 24, k)))
  theta_k <- Uk %*% diag(seq(0.5, 0.25, length.out = k), k, k) %*% t(Uk)
  Ak <- matrix(0, 24, 24)
  Ak[upper.tri(Ak)] <- theta_k[upper.tri(theta_k)] + rnorm(276, sd = 0.25)
  Ak <- Ak + t(Ak)
  fk <- fit_network(Ak, k, "gaussian", sigma = 0.25,
                    control = list(learning_rate = 0.02), U_star = Uk)
  check(fk$diagnostics$converged && fk$debiasing$status == "success",
        paste("K", k, "full Gaussian estimation and correction"))
  check(identical(dim(fk$Lambda_debiased), c(k, k)) &&
        length(fk$lambda_v_debiased) == k * (k + 1) / 2,
        paste("K", k, "full corrected output dimensions"))
  equal(crossprod(fk$U_hat) / 24, diag(k), paste("K", k, "full-fit normalization"))
}
Q <- matrix(c(cos(0.37), sin(0.37), -sin(0.37), cos(0.37)), 2)
Ur <- fit$U_hat %*% Q; Lr <- t(Q) %*% fit$Lambda_hat %*% Q
rotated <- api$.gfm_calculate_debiasing_components(Ur, Lr, U %*% Q, fit$config, A)
Lr_db <- api$.gfm_lambda_matrix(api$.gfm_get_lambda_v(Lr) + rotated$correction_total, K)
equal(Ur %*% Lr %*% t(Ur), fit$U_hat %*% fit$Lambda_hat %*% t(fit$U_hat),
      "raw signal is rotation invariant")
equal(Ur %*% Lr_db %*% t(Ur), fit$U_hat %*% fit$Lambda_debiased %*% t(fit$U_hat),
      "corrected signal is rotation invariant")

rotated_fit <- fit_network(A, K, "gaussian", sigma = 0.5, U_star = U %*% Q)
equal(rotated_fit$U_hat, Ur, "public API aligns to rotated truth")
equal(rotated_fit$Lambda_debiased, Lr_db, "public correction follows true coordinates")

# Missing, malformed, or mismatched truth must fail before fitting.
error(fit_network(A, K, "gaussian", 0.5), "requires the true U_star")
for (invalid_truth in list(U[, 1, drop = FALSE], U[-1, ], as.data.frame(U), U * NA)) {
  error(fit_network(A, K, "gaussian", 0.5, U_star = invalid_truth), "finite numeric n-by-K")
}
error(fit_network(A, K, "gaussian", 0.5, U_star = 2 * U), "crossprod")
U_named <- U; rownames(U_named) <- rev(rownames(A))
error(fit_network(A, K, "gaussian", 0.5, U_star = U_named), "same order")
rownames(U_named) <- rownames(A)
equal(api$.gfm_validate_U_star(U_named, A, K), U, "valid labeled truth retains values")
check(identical(rownames(api$.gfm_validate_U_star(U_named, A, K)), rownames(A)),
      "valid labeled truth retains row labels")

# Public input and control failures must occur before expensive fitting.
error(fit_network(as.data.frame(A), 2, "logit"), "numeric square")
error(fit_network(A[, -1], 2, "logit"), "numeric square")
B <- A; B[1, 2] <- NA; error(fit_network(B, 2, "gaussian", 1), "finite and observed")
B <- A; B[1, 2] <- B[1, 2] + 1; error(fit_network(B, 2, "gaussian", 1), "symmetric")
B <- A; colnames(B) <- rev(colnames(B)); error(fit_network(B, 2, "gaussian", 1), "same order")
error(fit_network(A, 2, "logit"), "binary")
error(fit_network(A, 2, "probit"), "binary")
error(fit_network(A, 2, "poisson"), "nonnegative integer")
B <- abs(A); error(fit_network(B, 2, "poisson"), "nonnegative integer")
B <- round(abs(A)); B[1, 2] <- B[2, 1] <- -1
error(fit_network(B, 2, "poisson"), "nonnegative integer")
error(fit_network(A, 2, "gaussian"), "Gaussian sigma")
error(fit_network(A, 2, "gaussian", sigma = 0), "Gaussian sigma")
error(fit_network(A, 2, "logit", sigma = 1), "only used")
error(fit_network(A, n, "gaussian", 1), "K must be smaller")
error(fit_network(A, 1.5, "gaussian", 1), "positive integer")
error(fit_network(A, 2, "gaussian", 1, debias = NA), "debias")
for (bad in list(list(unknown = 1), list(max_iter = 1.5), list(tolerance = 0),
                 list(probability_eps = 0.5), list(M_bound_U = 1), list(1),
                 list(learning_rate = 1, learning_rate = 2))) {
  error(fit_network(A, 2, "gaussian", 1, control = bad), "control|Unknown")
}
B <- A; diag(B) <- rep(c(NA, Inf, -999), length.out = n)
fit_diagonal <- fit_network(B, K, "gaussian", sigma = 0.5, debias = FALSE, U_star = U)
equal(fit_diagonal$Lambda_hat, fit$Lambda_hat, "diagonal ignored before initialization")
check(fit_diagonal$debiasing$status == "not_requested", "debias=FALSE status")
error(api$.gfm_project_C_U_certified(diag(5)[, 1:2], sqrt(2)), "joint C_U projection unavailable")

# Check CSV values and node order using temporary fixtures.
plain_csv <- tempfile(fileext = ".csv"); named_csv <- tempfile(fileext = ".csv")
write.table(A, plain_csv, sep = ",", row.names = FALSE, col.names = FALSE)
write.csv(A, named_csv)
equal(api$gfm_read_matrix(plain_csv), unname(A), "headerless CSV values")
named <- api$gfm_read_matrix(named_csv, node_names = TRUE)
equal(named, A, "named CSV values")
check(identical(dimnames(named), dimnames(A)), "named CSV node labels")
error(api$gfm_read_matrix(plain_csv, node_names = NA), "node_names")
B <- A; B[1, 2] <- B[1, 2] + 1
write.table(B, plain_csv, sep = ",", row.names = FALSE, col.names = FALSE)
error(api$gfm_read_matrix(plain_csv), "symmetric")
unlink(c(plain_csv, named_csv))

nonconverged <- warnings_from(fit_network(A, K, "gaussian", 0.5,
  control = list(max_iter = 1, tolerance = 1e-14), U_star = U))
check(length(nonconverged$messages) > 0 &&
      nonconverged$value$debiasing$status == "skipped", "nonconvergence retains raw fit and warns")
error(predict_network(nonconverged$value), "converged fit")
limited <- warnings_from(fit_network(A, K, "gaussian", 0.5,
  control = list(max_debias_memory_mb = 0.0001), U_star = U))
check(length(limited$messages) > 0 && limited$value$debiasing$status == "failed" &&
      grepl("workspace", limited$value$debiasing$message), "memory guard retains raw fit and warns")
equal(limited$value$Lambda_hat, fit$Lambda_hat, "memory guard leaves raw fit intact")
error(predict_network(limited$value, debiased = TRUE), "unavailable")
degenerate <- warnings_from(fit_network(matrix(0, 12, 12), 2, "gaussian", sigma = 1,
  U_star = sqrt(12) * qr.Q(qr(cbind(rep(1, 12), seq_len(12))))))
check(length(degenerate$messages) > 0 && degenerate$value$diagnostics$converged &&
      degenerate$value$debiasing$status == "failed" &&
      grepl("Remaining eigenvalues must be strictly positive",
            degenerate$value$debiasing$message),
      "zero network retains raw fit but rejects unidentified correction")
check(is.null(degenerate$value$Lambda_debiased), "rank failure has no corrected estimate")
prediction <- predict_network(fit, type = "link", debiased = TRUE)
expected <- fit$U_hat %*% fit$Lambda_debiased %*% t(fit$U_hat); diag(expected) <- NA_real_
equal(prediction, expected, "corrected prediction uses matching Lambda")
check(all(is.na(diag(prediction))) && identical(dimnames(prediction), dimnames(A)),
      "predictions exclude self-loops and preserve node order")
equal(predict_network(fit, type = "response"), predict_network(fit, type = "link"),
      "Gaussian response is its link")
error(predict_network(fit, debiased = NA), "debiased")
cat(sprintf("PASS: %d deterministic checks; no Monte Carlo simulations were run.\n", checks))
