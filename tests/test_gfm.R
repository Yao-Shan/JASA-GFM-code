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

# Parse only R code chunks, then evaluate ONLY selected function definitions.
# No library(), RNGkind(), truth generation, or runtime expressions are run.
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

# One observed-network fixture, not a Monte Carlo experiment.
n <- 35L; K <- 2L
U <- sqrt(n) * qr.Q(qr(matrix(rnorm(n * K), n, K)))
theta <- U %*% diag(c(0.5, -0.4)) %*% t(U)
A <- matrix(0, n, n)
A[upper.tri(A)] <- rnorm(n * (n - 1) / 2, theta[upper.tri(theta)], 0.5)
A <- A + t(A)
dimnames(A) <- rep(list(paste0("node_", seq_len(n))), 2)
before <- snapshot()
fit <- fit_network(A, K, "gaussian", sigma = 0.5)
check(identical(before, snapshot()), "fit has no RNG/options/search side effects")
check(fit$diagnostics$converged && fit$debiasing$status == "success", "Gaussian fit and correction succeed")
equal(crossprod(fit$U_hat) / n, diag(K), "fitted U normalization")
equal(fit$lambda_v_debiased, fit$lambda_v_hat + fit$correction$correction_total,
      "additive correction sign")
equal(api$.gfm_get_lambda_v(fit$Lambda_debiased), fit$lambda_v_debiased,
      "corrected matrix/vector correspondence")

# Initialization, finite-step optimization, and correction match original math
# for every likelihood. Probit thereby tests observed, outcome-dependent weights.
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
                                                  fit$U_hat, config, B)
  old <- reference$calculate_debiasing_components(fit$U_hat, fit$Lambda_hat,
                                                 fit$U_hat, config, B)
  equal(modern, old, paste(model, "correction parity in same fitted frame"))
}
raw <- reference$run_pgd_network(A, fit$config)
equal(raw$U_hat, fit$U_hat, "converged optimizer U parity")
equal(raw$Lambda_hat, fit$Lambda_hat, "converged optimizer Lambda parity")

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
                    control = list(learning_rate = 0.02))
  check(fk$diagnostics$converged && fk$debiasing$status == "success",
        paste("K", k, "full Gaussian estimation and correction"))
  check(identical(dim(fk$Lambda_debiased), c(k, k)) &&
        length(fk$lambda_v_debiased) == k * (k + 1) / 2,
        paste("K", k, "full corrected output dimensions"))
  equal(crossprod(fk$U_hat) / 24, diag(k), paste("K", k, "full-fit normalization"))
}
Q <- matrix(c(cos(0.37), sin(0.37), -sin(0.37), cos(0.37)), 2)
Ur <- fit$U_hat %*% Q; Lr <- t(Q) %*% fit$Lambda_hat %*% Q
rotated <- api$.gfm_calculate_debiasing_components(Ur, Lr, Ur, fit$config, A)
Lr_db <- api$.gfm_lambda_matrix(api$.gfm_get_lambda_v(Lr) + rotated$correction_total, K)
equal(Ur %*% Lr %*% t(Ur), fit$U_hat %*% fit$Lambda_hat %*% t(fit$U_hat),
      "raw signal is rotation invariant")
equal(Ur %*% Lr_db %*% t(Ur), fit$U_hat %*% fit$Lambda_debiased %*% t(fit$U_hat),
      "corrected signal is rotation invariant")

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
fit_diagonal <- fit_network(B, K, "gaussian", sigma = 0.5, debias = FALSE)
equal(fit_diagonal$Lambda_hat, fit$Lambda_hat, "diagonal ignored before initialization")
check(fit_diagonal$debiasing$status == "not_requested", "debias=FALSE status")
error(api$.gfm_project_C_U_certified(diag(5)[, 1:2], sqrt(2)), "joint C_U projection unavailable")
error(api$.gfm_safe_pinv_psd(diag(c(1, 0)), require_rank = TRUE), "unidentified directions")

# CSV imports preserve values and node order; no fixture is written to project.
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
  control = list(max_iter = 1, tolerance = 1e-14)))
check(length(nonconverged$messages) > 0 &&
      nonconverged$value$debiasing$status == "skipped", "nonconvergence retains raw fit and warns")
error(predict_network(nonconverged$value), "converged fit")
limited <- warnings_from(fit_network(A, K, "gaussian", 0.5,
  control = list(max_debias_memory_mb = 0.0001)))
check(length(limited$messages) > 0 && limited$value$debiasing$status == "failed" &&
      grepl("workspace", limited$value$debiasing$message), "memory guard retains raw fit and warns")
equal(limited$value$Lambda_hat, fit$Lambda_hat, "memory guard leaves raw fit intact")
error(predict_network(limited$value, debiased = TRUE), "unavailable")
degenerate <- warnings_from(fit_network(matrix(0, 12, 12), 2, "gaussian", sigma = 1))
check(length(degenerate$messages) > 0 && degenerate$value$diagnostics$converged &&
      degenerate$value$debiasing$status == "failed" &&
      grepl("information", degenerate$value$debiasing$message),
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
