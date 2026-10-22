# Validate saved-result semantics without simulations or a LaTeX installation.
args <- commandArgs(trailingOnly = FALSE)
this <- sub("^--file=", "", args[grepl("^--file=", args)])
this <- gsub("~+~", " ", this, fixed = TRUE)
root <- normalizePath(file.path(dirname(this), ".."), mustWork = TRUE)
source(file.path(root, "plots", "plot_individual.R"))
stop_if_error_missing <- function(expr) stopifnot(inherits(tryCatch({force(expr); NULL}, error = identity), "error"))
work <- tempfile("gfm_plotting_test_"); dir.create(work)
for (K in c(1L, 3L)) {
  folder <- file.path(work, paste0("K", K)); dir.create(folder)
  P <- K * (K + 1L) / 2L
  cfg <- list(model = "gaussian", n = 10L, k = K, n_sim = 4L,
    lambda_vec = seq_len(K) / 5, target_node = 1L, design_seed = 1L,
    replication_seed_base = 10L, truth_design = "test_fixture")
  u <- matrix(seq_len(4L * K) / 20, 4L, K)
  z <- matrix(seq_len(4L * P) / 10, 4L, P)
  e <- matrix(seq_len(4L * P) / 100, 4L, P)
  d <- data.frame(convergence_code = rep(0L, 4L), chi_sq_U = rowSums(u^2),
                  chi_sq_debiased_lambda = rowSums(z^2))
  for (j in seq_len(K)) d[[paste0("stat_U_", j)]] <- u[, j]
  for (j in seq_len(P)) {
    d[[paste0("debiased_lambda_", j)]] <- z[, j]
    d[[paste0("lambda_raw_error_", j)]] <- e[, j]
    d[[paste0("lambda_debiased_error_", j)]] <- e[, j] - .01
  }
  csv <- file.path(folder, "case_results.csv"); config <- file.path(folder, "case_config.rds")
  write.csv(d, csv, row.names = FALSE); saveRDS(cfg, config)
  x <- gfm_read_figure_cases(folder, n = 10, models = "gaussian")[[1]]
  stopifnot(x$K == K, x$P == P,
    isTRUE(all.equal(unname(x$raw_error[, seq_len(K), drop = FALSE]), e[, seq_len(K), drop = FALSE] * sqrt(2))),
    isTRUE(all.equal(unname(x$raw_estimate - x$raw_error), matrix(x$true_entries, 4L, P, byrow = TRUE))),
    identical(unname(as.matrix(x$d[paste0("debiased_lambda_", seq_len(P))])), unname(as.matrix(d[paste0("debiased_lambda_", seq_len(P))]))))
  if (K == 3L) stopifnot(identical(as.integer(x$coordinates$a), c(1L,2L,3L,1L,1L,2L)),
                          identical(as.integer(x$coordinates$b), c(1L,2L,3L,2L,3L,3L)))
  bad <- d; bad$lambda_raw_error_1 <- NULL; write.csv(bad, csv, row.names = FALSE)
  stop_if_error_missing(gfm_read_figure_cases(folder, n = 10))
  stopifnot(length(gfm_read_figure_cases(folder, n = 10, figure_types = "qq")) == 1L)
  histogram_only <- d[!grepl("^(stat_U_|debiased_lambda_|chi_sq_)", names(d))]
  write.csv(histogram_only, csv, row.names = FALSE)
  stopifnot(length(gfm_read_figure_cases(folder, n = 10, figure_types = "hist_errors")) == 1L)
  bad <- d; bad$convergence_code[1] <- 1L; write.csv(bad, csv, row.names = FALSE)
  stop_if_error_missing(gfm_read_figure_cases(folder, n = 10))
  retained <- suppressWarnings(gfm_read_figure_cases(folder, n = 10, require_complete = FALSE)[[1]])
  stopifnot(nrow(retained$d) == 3L, sum(!retained$selection$included) == 1L)
  write.csv(d, csv, row.names = FALSE)
  file.copy(csv, file.path(folder, "duplicate_results.csv"))
  file.copy(config, file.path(folder, "duplicate_config.rds"))
  stop_if_error_missing(gfm_read_figure_cases(folder, n = 10))
}
unlink(work, recursive = TRUE)
cat("PASS: K=1/K=3 Lambda units and ordering; unmodified standardized QQ statistics; type-specific required columns; recorded exclusions; missing errors, failed rows, and duplicate runs rejected.\n")
