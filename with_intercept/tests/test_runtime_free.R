# Compare ordinary-R sequential and PSOCK results on the same two networks.
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
stopifnot(length(arg) == 1L)
root <- dirname(dirname(normalizePath(gsub("~+~", " ",
  sub("^--file=", "", arg), fixed = TRUE), mustWork = TRUE)))
stopifnot(!dir.exists(file.path(root, "runtime")))
api <- new.env(parent = globalenv())
source(file.path(root, "R", "gfm_baseline.R"), local = api)
arguments <- list(n = 24L, K = 2L, models = "gaussian", n_sim = 2L,
  nodes = 1L, alpha_seed = 817285L, replication_seed_base = 2000L,
  control = list(tolerance = 0.05), debias = TRUE,
  inference_workers = 1L, output_dir = NULL, verbose = FALSE)
sequential <- do.call(api$gfb_run_simulation,
  c(arguments, list(n_cores = 1L, worker_backend = "in_process")))
parallel <- do.call(api$gfb_run_simulation,
  c(arguments, list(n_cores = 2L, worker_backend = "psock")))
for (result in list(sequential, parallel)) {
  stopifnot(length(result$records) == 2L,
    !result$settings$use_fast_fit, !result$settings$use_fast_inference)
  for (backend in result$settings$runtime_preflight) {
    stopifnot(identical(backend$status, "disabled"),
              identical(backend$backend, "pure_R"))
  }
  for (worker in result$settings$worker_runtime) for (backend in worker$runtime) {
    stopifnot(identical(backend$status, "disabled"),
              identical(backend$backend, "pure_R"))
  }
  stopifnot(all(result$fit_summary$status == "success"),
            all(result$lambda_summary$status == "ok"))
}
stopifnot(isTRUE(all.equal(sequential$node_summary, parallel$node_summary,
                          tolerance = 1e-8)),
          isTRUE(all.equal(sequential$lambda_summary, parallel$lambda_summary,
                          tolerance = 1e-8)))
cat("PASS: sequential and two-worker PSOCK fitting/inference used pure R without runtime files; results agree.\n")
