# Monte Carlo driver for the GFM with additive node effects.
# Load through gfm_baseline.R; sourcing alone does not run a study.

.gfb_bind_rows <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(data.frame())
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(x) {
    for (nm in setdiff(columns, names(x))) x[[nm]] <- NA
    x[, columns, drop = FALSE]
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

.gfb_sim_scalar <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) != 1L) default else x
}

.gfb_run_directory <- function(output_dir, run_id) {
  if (is.null(output_dir)) return(NULL)
  if (!is.character(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) stop("output_dir must be a directory path or NULL.")
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Could not create output_dir: ", output_dir)
  }
  if (!is.null(run_id)) {
    if (!is.character(run_id) || length(run_id) != 1L || is.na(run_id) ||
        !grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", run_id) || run_id %in% c(".", "..")) {
      stop("run_id must be one folder name containing letters, numbers, '.', '_' or '-'.")
    }
    candidate <- file.path(output_dir, run_id)
    if (file.exists(candidate)) stop("Run already exists; refusing to overwrite: ", candidate)
    if (!dir.create(candidate)) stop("Could not create run directory: ", candidate)
  } else {
    stem <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
    suffix <- 0L
    repeat {
      candidate <- file.path(output_dir, paste0(stem, if (suffix) paste0("_", suffix) else ""))
      if (!file.exists(candidate) && dir.create(candidate, showWarnings = FALSE)) break
      suffix <- suffix + 1L
      if (suffix > 10000L) stop("Could not reserve a new run directory.")
    }
  }
  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

.gfb_sim_summaries <- function(records, K) {
  fit_rows <- node_rows <- lambda_rows <- vector("list", length(records))
  for (ii in seq_along(records)) {
    rec <- records[[ii]]
    if (is.null(rec)) next
    diag <- rec$diagnostics
    inf <- rec$inference
    fit_rows[[ii]] <- data.frame(
      model = rec$model, replicate = rec$replicate, network_seed = rec$network_seed,
      status = rec$status, converged = isTRUE(rec$converged),
      convergence_code = .gfb_sim_scalar(diag$convergence_code),
      n_iter = .gfb_sim_scalar(diag$n_iter),
      scaled_gradient = .gfb_sim_scalar(diag$scaled_gradient),
      center_error = .gfb_sim_scalar(diag$center_error),
      gram_error = .gfb_sim_scalar(diag$gram_error),
      max_row_norm = .gfb_sim_scalar(diag$max_row_norm),
      alpha_bound_margin = .gfb_sim_scalar(diag$alpha_bound_margin),
      one_sided_node_count = length(diag$one_sided_nodes),
      inference_status = .gfb_sim_scalar(inf$status, NA_character_),
      lambda_status = .gfb_sim_scalar(inf$lambda$status, NA_character_),
      elapsed_seconds = rec$elapsed_seconds,
      fit_elapsed_seconds = .gfb_sim_scalar(rec$fit_elapsed_seconds),
      inference_elapsed_seconds = .gfb_sim_scalar(rec$inference_elapsed_seconds),
      inference_wait_seconds = .gfb_sim_scalar(rec$inference_wait_seconds),
      worker_pid = .gfb_sim_scalar(rec$worker$pid),
      error = .gfb_sim_scalar(rec$error, NA_character_),
      warnings = paste(rec$warnings, collapse = " | "), stringsAsFactors = FALSE
    )
    if (is.data.frame(inf$nodes) && nrow(inf$nodes)) {
      node_rows[[ii]] <- cbind(
        data.frame(model = rec$model, replicate = rec$replicate,
                   fit_converged = isTRUE(rec$converged), stringsAsFactors = FALSE),
        inf$nodes
      )
    }
    lam <- inf$lambda
    if (!is.null(lam)) {
      m <- K * (K + 1L) / 2L
      row <- data.frame(
        model = rec$model, replicate = rec$replicate,
        fit_converged = isTRUE(rec$converged),
        status = .gfb_sim_scalar(lam$status, NA_character_),
        debias_method = .gfb_sim_scalar(lam$method,
          .gfb_sim_scalar(rec$debias_method, NA_character_)),
        message = .gfb_sim_scalar(lam$message, ""),
        chi_sq_raw = .gfb_sim_scalar(lam$chi_sq_raw),
        chi_sq_debiased = .gfb_sim_scalar(lam$chi_sq_debiased),
        df = .gfb_sim_scalar(lam$df, m), stringsAsFactors = FALSE
      )
      fields <- c("lambda_v_hat", "lambda_v_debiased", "raw_error", "debiased_error",
                  "z_raw", "z_debiased", "correction_total",
                  "correction_likelihood", "correction_ell3", "correction_score_ell2", "correction_BM",
                  "correction_constraint")
      for (field in fields) {
        values <- lam[[field]]
        if (is.null(values) && field == "correction_BM") values <- lam$correction_score_ell2
        if (is.null(values) && field == "correction_score_ell2") values <- lam$correction_BM
        if (is.null(values) || length(values) != m) values <- rep(NA_real_, m)
        for (jj in seq_len(m)) row[[paste0(field, "_", jj)]] <- values[jj]
      }
      lambda_rows[[ii]] <- row
    }
  }
  list(fit_summary = .gfb_bind_rows(fit_rows),
       node_summary = .gfb_bind_rows(node_rows),
       lambda_summary = .gfb_bind_rows(lambda_rows))
}

# Export named R functions; install optional native functions inside each worker.
.gfb_worker_function_names <- function() {
  c(".gfb_normalize_model_name", ".gfb_symmetrize", ".gfb_sigmoid",
    ".gfb_stable_log1pexp", ".gfb_safe_pinv_psd", ".gfb_assert_positive_definite",
    ".gfb_ordinary_inverse_pd", ".gfb_get_lambda_v", ".gfb_lambda_direction_derivative",
    ".gfb_proximal_projection_Lambda", ".gfb_probit_ratios", ".gfb_link_probability",
    ".gfb_link_inverse", ".gfb_negative_loglikelihood_edges", ".gfb_negative_score_matrix",
    ".gfb_expected_derivative_terms", ".gfb_observed_information_matrix",
    ".gfb_spectral_keep_indices", ".gfb_validate_matrix", ".gfb_lambda_matrix",
    ".gfb_scalar", "gfb_config", ".gfb_with_seed", ".gfb_theta", ".gfb_safe_sqrt_psd",
    ".gfb_prepare_alpha",
    "gfb_generate_truth", "gfb_generate_network", ".gfb_project_U", ".gfb_objective",
    ".gfb_gradient", ".gfb_additive_initial", ".gfb_initialize", ".gfb_step", "gfb_fit",
    "gfb_predict", "print.gfb_fit", "gfb_align", ".gfb_build_D_YG", ".gfb_debias_memory_mb",
    ".gfb_information", ".gfb_tangent_covariance", ".gfb_w_lambda_derivative",
    ".gfb_bias_terms", ".gfb_bias_terms_pdf_bp_bm", "gfb_infer", ".gfb_atomic_rds", ".gfb_sim_progress",
    ".gfb_inference_slot", ".gfb_sim_task", ".gfb_worker_run", ".gfb_runtime_install")
}

.gfb_atomic_rds <- function(value, path, overwrite = FALSE) {
  tmp <- tempfile(pattern = ".gfb-writing-", tmpdir = dirname(path))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(value, tmp)
  if (overwrite) {
    if (!file.rename(tmp, path)) stop("Could not atomically replace ", path)
  } else {
    if (!file.link(tmp, path)) stop("Refusing to overwrite or unable to create ", path)
  }
  invisible(path)
}

.gfb_sim_progress <- function(task, phase, directory, extra = list()) {
  if (is.null(directory)) return(invisible(NULL))
  .gfb_atomic_rds(c(list(index = task$index, model = task$model,
    replicate = task$replicate, network_seed = task$network_seed, phase = phase,
    pid = Sys.getpid(), updated_at = as.character(Sys.time())), extra),
    file.path(directory, sprintf("progress_%05d.rds", task$index)), overwrite = TRUE)
}

# Reserve an inference slot atomically; only its owner releases it.
# Interrupted workers may leave stale slots.
.gfb_inference_slot <- function(directory, n_slots, task, progress_dir, expr) {
  started <- proc.time()[["elapsed"]]
  slot <- NULL
  repeat {
    for (i in seq_len(n_slots)) {
      candidate <- file.path(directory, paste0("slot_", i))
      if (dir.create(candidate, showWarnings = FALSE)) {
        slot <- candidate
        break
      }
    }
    if (!is.null(slot)) break
    Sys.sleep(.1)
  }
  owner <- list(pid = Sys.getpid(), index = task$index,
                acquired_at = as.character(Sys.time()))
  on.exit({
    saved <- tryCatch(readRDS(file.path(slot, "owner.rds")), error = function(e) NULL)
    if (identical(saved, owner)) unlink(slot, recursive = TRUE)
  }, add = TRUE)
  .gfb_atomic_rds(owner, file.path(slot, "owner.rds"))
  wait <- proc.time()[["elapsed"]] - started
  .gfb_sim_progress(task, "inference", progress_dir,
                    list(inference_wait_seconds = wait, inference_slot = basename(slot)))
  value <- force(expr)
  list(value = value, wait_seconds = unname(wait))
}

.gfb_runtime_install <- function(envir, runtime_dir, use_fast_fit,
                                 use_fast_inference, debias_method) {
  result <- list(fit = list(status = "disabled", backend = "pure_R"),
                 inference = list(status = "disabled", backend = "pure_R"))
  # Serialize optional native installation to avoid cache races; fits stay parallel.
  lock_wait_seconds <- 0
  needs_runtime <- isTRUE(use_fast_fit) ||
    (isTRUE(use_fast_inference) && debias_method %in% c("pdf_bp_bm", "constrained_second_order"))
  if (needs_runtime && !is.null(runtime_dir) && dir.exists(runtime_dir)) {
    lock <- file.path(runtime_dir, ".gfb-runtime-install-lock")
    owner <- list(pid = Sys.getpid(), acquired_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6%z"),
      process_elapsed = unname(proc.time()[["elapsed"]]))
    wait_started <- proc.time()[["elapsed"]]
    repeat {
      if (dir.create(lock, showWarnings = FALSE)) break
      if (proc.time()[["elapsed"]] - wait_started >= 120)
        stop("Timed out waiting for the runtime installation lock: ", lock,
          ". An active or unknown lock is never removed automatically; inspect owner.rds before retrying.")
      Sys.sleep(.1)
    }
    on.exit({
      saved <- tryCatch(readRDS(file.path(lock, "owner.rds")), error = function(e) NULL)
      if (identical(saved, owner)) unlink(lock, recursive = TRUE)
    }, add = TRUE)
    .gfb_atomic_rds(owner, file.path(lock, "owner.rds"))
    lock_wait_seconds <- unname(proc.time()[["elapsed"]] - wait_started)
  }
  install <- function(kind, installer) {
    path <- if (is.null(runtime_dir)) "" else file.path(runtime_dir, paste0("fast_", kind, ".R"))
    if (!nzchar(path) || !file.exists(path)) return(list(status = "fallback",
      backend = "pure_R", reason = paste("Missing fast", kind, "runtime file")))
    original <- lapply(c(".gfb_objective", ".gfb_gradient", ".gfb_information", ".gfb_bias_terms"),
      function(nm) get(nm, envir = envir, inherits = FALSE))
    names(original) <- c(".gfb_objective", ".gfb_gradient", ".gfb_information", ".gfb_bias_terms")
    tryCatch({
      runtime <- new.env(parent = globalenv())
      source(path, local = runtime)
      installed <- installer(runtime)
      result <- list(status = "installed", backend = paste0("Rcpp_fast_", kind))
      if (is.list(installed) && !is.null(installed$reported_backend)) {
        result$backend <- installed$reported_backend
        result$fused_status <- installed$fused_status
        if (!is.null(installed$reason)) result$reason <- installed$reason
      }
      result
    }, error = function(e) {
      list2env(original, envir = envir)
      list(status = "fallback", backend = "pure_R", reason = conditionMessage(e))
    })
  }
  if (use_fast_fit) result$fit <- install("fit", function(runtime) {
    runtime$gfb_install_fast_fit(envir, cache_dir = file.path(runtime_dir, "cache", "fit"))
    fused_path <- file.path(runtime_dir, "fast_fit_fused.R")
    if (!file.exists(fused_path)) return(list(reported_backend = "Rcpp_fast_fit",
      fused_status = "not_available"))
    # Load the fused kernel from the explicitly selected runtime directory.
    fast_helpers <- mget(c(".gfb_objective", ".gfb_gradient"), envir = envir)
    fused <- tryCatch({
      source(fused_path, local = runtime)
      runtime$gfb_install_fused_fit(envir,
        cache_dir = file.path(runtime_dir, "cache", "fused_fit"))
      TRUE
    }, error = identity)
    if (inherits(fused, "error")) {
      list2env(fast_helpers, envir = envir)
      return(list(reported_backend = "Rcpp_fast_fit", fused_status = "fallback",
        reason = paste("Fused acceleration unavailable; retaining fast_fit:", conditionMessage(fused))))
    }
    list(reported_backend = "Rcpp_fused_fit", fused_status = "installed")
  })
  if (use_fast_inference && debias_method %in% c("pdf_bp_bm", "constrained_second_order")) {
    # Share exact information assembly; use this bias helper only for
    # constrained_second_order inference.
    result$inference <- install("inference", function(runtime) {
      runtime$gfb_install_fast_inference(envir,
        cpp_file = file.path(runtime_dir, "fast_inference.cpp"),
        cache_dir = file.path(runtime_dir, "cache", "inference"))
    })
  } else if (use_fast_inference) {
    result$inference <- list(status = "disabled_for_method", backend = "pure_R",
      reason = paste("Fast inference does not support method:", debias_method))
  }
  for (kind in c("fit", "inference")) result[[kind]]$installation_lock_wait_seconds <- lock_wait_seconds
  result
}

# Install truth, controls and functions once per worker in a private environment.
.gfb_prepare_worker <- function(fns, task_args, runtime_dir, use_fast_fit,
                                use_fast_inference, debias_method, paths,
                                inference_workers, envir = .GlobalEnv) {
  for (nm in names(fns)) {
    environment(fns[[nm]]) <- envir
    assign(nm, fns[[nm]], envir = envir)
  }
  assign(".gfb_root", if (is.null(runtime_dir)) getwd() else dirname(runtime_dir), envir)
  assign(".gfb_debias_method", debias_method, envir)
  Sys.setenv(VECLIB_MAXIMUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1")
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L)
    RhpcBLASctl::omp_set_num_threads(1L)
  }
  runtime <- get(".gfb_runtime_install", envir)(envir, runtime_dir,
    use_fast_fit, use_fast_inference, debias_method)
  info <- list(pid = Sys.getpid(), BLAS = unname(extSoftVersion()["BLAS"]),
    LAPACK = La_library(), threads = Sys.getenv(c("VECLIB_MAXIMUM_THREADS",
      "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")),
    runtime = runtime, debias_method = debias_method)
  assign(".gfb_worker_state", list(args = task_args, info = info, paths = paths,
    inference_workers = inference_workers), envir)
  info
}

.gfb_worker_run <- function(task) {
  state <- .gfb_worker_state
  rec <- do.call(.gfb_sim_task, c(list(task = task), state$args,
    list(worker = state$info, paths = state$paths, inference_workers = state$inference_workers)))
  if (!is.null(state$paths$partial)) .gfb_atomic_rds(rec,
    file.path(state$paths$partial, sprintf("record_%05d.rds", task$index)))
  rec
}

.gfb_sim_task <- function(task, truth, K, sigma, control, nodes, debias, save_fits,
                          worker = NULL, paths = list(progress = NULL, gate = NULL),
                          inference_workers = 1L,
                          debias_method = "pdf_bp_bm") {
  started <- proc.time()[["elapsed"]]
  rec <- list(model = task$model, replicate = task$replicate, network_seed = task$network_seed,
    status = "started", converged = FALSE, diagnostics = NULL, estimates = NULL,
    inference = NULL, warnings = character(), error = NULL, worker = worker,
    debias_method = if (debias) debias_method else "none",
    effective_control = gfb_config(nrow(truth$U_star), K, task$model, sigma, control),
    truth_seeds = truth$seeds, fit_elapsed_seconds = NA_real_,
    inference_elapsed_seconds = NA_real_, inference_wait_seconds = 0)
  phase <- "generation"
  phase_started <- started
  captured_warnings <- character()
  result <- tryCatch(withCallingHandlers({
    .gfb_sim_progress(task, phase, paths$progress)
    A <- gfb_generate_network(truth, model = task$model, sigma = sigma, seed = task$network_seed)
    phase <- "fit"; phase_started <- proc.time()[["elapsed"]]
    .gfb_sim_progress(task, phase, paths$progress)
    fit <- gfb_fit(A, K = K, model = task$model, sigma = sigma, control = control)
    rec$fit_elapsed_seconds <- unname(proc.time()[["elapsed"]] - phase_started)
    rec$diagnostics <- fit$diagnostics
    rec$converged <- isTRUE(fit$diagnostics$convergence_code == 0L)
    rec$estimates <- list(U_hat = fit$U_hat, alpha_hat = fit$alpha_hat, Lambda_hat = fit$Lambda_hat)
    if (save_fits) rec$fit <- fit
    if (!rec$converged) {
      rec$status <- "fit_not_converged"
      rec$error <- "Optimizer did not converge; inference was not computed."
    } else {
      phase <- "inference"; phase_started <- proc.time()[["elapsed"]]
      infer_once <- function() gfb_infer(fit, U_star = truth$U_star,
        alpha_star = truth$alpha_star, Lambda_star = truth$Lambda_star,
        nodes = nodes, debias = debias, debias_method = debias_method)
      .gfb_sim_progress(task, if (debias) "waiting_for_inference" else phase,
        paths$progress, list(fit_elapsed_seconds = rec$fit_elapsed_seconds,
          n_iter = fit$diagnostics$n_iter, scaled_gradient = fit$diagnostics$scaled_gradient))
      if (debias && !is.null(paths$gate)) {
        gated <- .gfb_inference_slot(paths$gate, inference_workers, task, paths$progress, infer_once())
        rec$inference <- gated$value
        rec$inference_wait_seconds <- gated$wait_seconds
      } else rec$inference <- infer_once()
      rec$inference_elapsed_seconds <- unname(proc.time()[["elapsed"]] - phase_started - rec$inference_wait_seconds)
      rec$status <- switch(rec$inference$status,
        ok = "success", failed = "inference_failed", "inference_partial")
      if (!identical(rec$inference$status, "ok")) rec$error <- rec$inference$message
    }
    NULL
  }, warning = function(w) {
    captured_warnings <<- c(captured_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(result, "error")) {
    rec$status <- paste0(phase, "_error")
    rec$error <- conditionMessage(result)
    if (phase == "fit") rec$fit_elapsed_seconds <- unname(proc.time()[["elapsed"]] - phase_started)
    if (phase == "inference") rec$inference_elapsed_seconds <- unname(proc.time()[["elapsed"]] - phase_started)
  }
  rec$warnings <- unique(captured_warnings)
  rec$elapsed_seconds <- unname(proc.time()[["elapsed"]] - started)
  .gfb_sim_progress(task, "finished", paths$progress,
    list(status = rec$status, fit_elapsed_seconds = rec$fit_elapsed_seconds,
      inference_elapsed_seconds = rec$inference_elapsed_seconds,
      inference_wait_seconds = rec$inference_wait_seconds, elapsed_seconds = rec$elapsed_seconds,
      diagnostics = rec$diagnostics, error = rec$error))
  rec
}

# Run independently seeded networks conditional on fixed truth.
# Preserve failed attempts and save each completed task in partial/.
gfb_run_simulation <- function(
    n = 60L, K = 2L, models = c("gaussian", "logit", "poisson", "probit"),
    n_sim = 3L, nodes = 1L, lambda_values = c(0.5, 0.5), sigma = 1,
    design_seed = 238687L, alpha_seed = 817263L,
    replication_seed_base = 2000L, alpha_min = -2, alpha_max = 0,
    control = list(), debias = TRUE, output_dir = NULL, run_id = NULL,
    save_fits = FALSE, verbose = TRUE, n_cores = 7L, use_fast_fit = FALSE,
    runtime_dir = NULL, inference_workers = 2L, use_fast_inference = FALSE,
    log_file = NULL, worker_backend = "auto",
    debias_method = c("pdf_bp_bm", "constrained_second_order"),
    alpha_normalization = c("l2", "none")) {
  debias_method <- match.arg(debias_method)
  alpha_normalization <- match.arg(alpha_normalization)
  integer_scalar <- function(x, name, minimum = 1L, maximum = .Machine$integer.max) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != floor(x) || x < minimum || x > maximum) stop(name, " must be an integer in ", minimum, "..", maximum, ".")
    as.integer(x)
  }
  n <- integer_scalar(n, "n", 3L); K <- integer_scalar(K, "K")
  n_sim <- integer_scalar(n_sim, "n_sim")
  n_cores <- integer_scalar(n_cores, "n_cores", maximum = 7L)
  inference_workers <- integer_scalar(inference_workers, "inference_workers", maximum = 2L)
  worker_backend <- match.arg(worker_backend, c("auto", "psock", "in_process"))
  replication_seed_base <- integer_scalar(replication_seed_base, "replication_seed_base", 0L)
  if (as.double(replication_seed_base) + n_sim > .Machine$integer.max)
    stop("replication_seed_base + n_sim exceeds R's integer seed range.")
  if (!is.numeric(nodes) || !length(nodes) || anyNA(nodes) || any(!is.finite(nodes)) ||
      any(nodes != floor(nodes)) || any(nodes < 1L | nodes > n) || anyDuplicated(nodes))
    stop("nodes must contain distinct node indices in 1,...,n.")
  nodes <- as.integer(nodes)
  for (flag in c("debias", "save_fits", "verbose", "use_fast_fit", "use_fast_inference")) {
    value <- get(flag)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) stop(flag, " must be TRUE or FALSE.")
  }
  if (!is.character(models) || !length(models) || anyNA(models)) stop("models must be a nonempty character vector.")
  configs <- lapply(models, function(model) gfb_config(n, K, model, sigma, control))
  models <- vapply(configs, function(config) config$model, character(1))
  if (anyDuplicated(models)) stop("models contains duplicate observation models after normalization.")
  names(configs) <- models
  n_workers <- min(n_cores, length(models) * n_sim)
  if (worker_backend == "in_process" && n_workers > 1L)
    stop("worker_backend='in_process' requires one worker; use n_cores=1.")
  inference_workers <- min(inference_workers, n_workers)
  api <- environment(gfb_run_simulation)
  names_to_copy <- .gfb_worker_function_names()
  names_to_copy <- names_to_copy[vapply(names_to_copy, exists, logical(1), envir = api, inherits = FALSE)]
  fns <- mget(names_to_copy, envir = api, inherits = FALSE)
  if (!all(vapply(fns, is.function, logical(1)))) stop("The worker export list contains a non-function.")
  # Reject native closures that cannot be transferred to fresh PSOCK workers.
  bad <- names(fns)[!vapply(fns, function(f) identical(environment(f), api), logical(1))]
  if (length(bad)) stop("Source a clean gfm_baseline.R before running; nonportable closures found: ", paste(bad, collapse = ", "))
  for (nm in names(fns)) environment(fns[[nm]]) <- globalenv()
  if (is.null(runtime_dir)) {
    root <- get0(".gfb_root", envir = api, inherits = FALSE)
    if (!is.null(root)) runtime_dir <- file.path(root, "runtime")
  }
  if (!is.null(runtime_dir)) {
    if (!is.character(runtime_dir) || length(runtime_dir) != 1L || is.na(runtime_dir)) stop("runtime_dir must be one directory or NULL.")
    runtime_dir <- normalizePath(runtime_dir, winslash = "/", mustWork = FALSE)
  }
  launcher <- if (is.null(runtime_dir)) "" else file.path(runtime_dir, "rscript_accelerate.sh")
  blas_dir <- if (is.null(runtime_dir)) "" else file.path(runtime_dir, "blas")
  use_accelerate_launcher <- identical(unname(Sys.info()["sysname"]), "Darwin") &&
    nzchar(launcher) && file.exists(launcher) && file.access(launcher, mode = 1L) == 0L &&
    file.exists(file.path(blas_dir, "libRblas.dylib"))
  # Use a fresh R process for PSOCK, including single-worker runs.
  actual_worker_backend <- if (worker_backend == "auto") {
    if (n_workers > 1L || use_accelerate_launcher) "psock" else "in_process"
  } else worker_backend
  # Warm optional caches without exporting compiled closures.
  precompile <- new.env(parent = globalenv())
  for (nm in names(fns)) {
    fn <- fns[[nm]]; environment(fn) <- precompile; assign(nm, fn, precompile)
  }
  assign(".gfb_root", if (is.null(runtime_dir)) getwd() else dirname(runtime_dir), precompile)
  runtime_preflight <- .gfb_runtime_install(precompile, runtime_dir, use_fast_fit,
    use_fast_inference && debias, debias_method)
  # Fall back once if parent compilation fails.
  worker_fast_fit <- use_fast_fit && identical(runtime_preflight$fit$status, "installed")
  worker_fast_inference <- use_fast_inference && debias && identical(runtime_preflight$inference$status, "installed")
  truth <- gfb_generate_truth(n, K, lambda_values, design_seed, alpha_seed, alpha_min, alpha_max,
                            alpha_normalization = alpha_normalization)
  settings <- list(schema_version = "gfm-baseline-v1", n = n, K = K, n_sim = n_sim,
    models = models, nodes = nodes, lambda_values = lambda_values, sigma = sigma,
    design_seed = design_seed, alpha_seed = alpha_seed, replication_seed_base = replication_seed_base,
    network_seeds = replication_seed_base + seq_len(n_sim), alpha_min = alpha_min, alpha_max = alpha_max,
    alpha_normalization = alpha_normalization,
    alpha_raw_distribution = truth$alpha_raw_distribution,
    alpha_distribution = truth$alpha_distribution, alpha_raw_norm = truth$alpha_raw_norm,
    alpha_normalization_metadata = truth$alpha_normalization_metadata,
    control = control, configs = configs, debias = debias,
    debias_method = if (debias) debias_method else "none", save_fits = save_fits,
    execution = if (n_workers > 1L) "parallel" else "sequential", n_workers = n_workers,
    worker_backend_requested = worker_backend, worker_backend = actual_worker_backend,
    inference_workers = inference_workers, fixed_truth = TRUE,
    use_fast_fit = use_fast_fit, use_fast_inference = use_fast_inference,
    runtime_dir = runtime_dir, runtime_preflight = runtime_preflight,
    seed_policy = "replication_seed_base + replicate; reused separately across models",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), session_info = utils::sessionInfo())
  run_dir <- .gfb_run_directory(output_dir, run_id)
  scratch <- if (is.null(run_dir)) tempfile("gfb-run-") else run_dir
  if (is.null(run_dir)) dir.create(scratch)
  paths <- list(partial = if (is.null(run_dir)) NULL else file.path(run_dir, "partial"),
    progress = file.path(scratch, "progress"), gate = file.path(scratch, "inference_gate"))
  for (path in Filter(Negate(is.null), paths)) dir.create(path)
  if (!is.null(run_dir)) {
    .gfb_atomic_rds(settings, file.path(run_dir, "settings.rds"))
    .gfb_atomic_rds(truth, file.path(run_dir, "truth.rds"))
  }
  if (is.null(log_file) && !is.null(run_dir)) log_file <- file.path(run_dir, "progress.log")
  if (!is.null(log_file) && (!is.character(log_file) || length(log_file) != 1L || is.na(log_file))) stop("log_file must be one path or NULL.")
  log_line <- function(line) {
    line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), line)
    if (verbose) { message(line); flush.console() }
    if (!is.null(log_file)) cat(line, "\n", file = log_file, append = TRUE)
  }
  log_line(sprintf("%d tasks; %d workers; backend=%s; at most %d concurrent dense inferences; method=%s",
    length(models) * n_sim, n_workers, actual_worker_backend, inference_workers, settings$debias_method))
  for (model in models) log_line(sprintf("%s: learning_rate=%g; tolerance=%g; max_iter=%d",
    model, configs[[model]]$learning_rate, configs[[model]]$tolerance, configs[[model]]$max_iter))
  for (kind in names(runtime_preflight)) log_line(paste("runtime", kind,
    runtime_preflight[[kind]]$status, .gfb_sim_scalar(runtime_preflight[[kind]]$reason, "")))
  tasks <- vector("list", length(models) * n_sim); rr <- 0L
  for (model in models) for (replicate in seq_len(n_sim)) {
    rr <- rr + 1L
    tasks[[rr]] <- list(index = rr, model = model, replicate = replicate,
      network_seed = replication_seed_base + replicate)
  }
  task_args <- list(truth = truth, K = K, sigma = sigma, control = control,
    nodes = nodes, debias = debias, save_fits = save_fits, debias_method = debias_method)
  records <- vector("list", length(tasks)); completed <- 0L
  study_started <- proc.time()[["elapsed"]]
  record_result <- function(rec, index) {
    if (inherits(rec, "try-error") || !is.list(rec) || is.null(rec$status)) {
      rec <- list(model = tasks[[index]]$model, replicate = tasks[[index]]$replicate,
        network_seed = tasks[[index]]$network_seed, status = "worker_error", converged = FALSE,
        diagnostics = NULL, estimates = NULL, inference = NULL, warnings = character(),
        error = paste(rec, collapse = " "), elapsed_seconds = NA_real_)
      if (!is.null(paths$partial)) .gfb_atomic_rds(rec,
        file.path(paths$partial, sprintf("record_%05d.rds", index)))
    }
    records[[index]] <<- rec; completed <<- completed + 1L
    log_line(sprintf("[%s] replicate %d/%d: %s; iter=%s; residual=%s; fit=%.2fs; infer=%.2fs; gate=%.2fs; task=%.2fs; completed=%d/%d; wall=%.2fs",
      rec$model, rec$replicate, n_sim, rec$status,
      .gfb_sim_scalar(rec$diagnostics$n_iter, NA), .gfb_sim_scalar(rec$diagnostics$scaled_gradient, NA),
      .gfb_sim_scalar(rec$fit_elapsed_seconds), .gfb_sim_scalar(rec$inference_elapsed_seconds),
      .gfb_sim_scalar(rec$inference_wait_seconds), rec$elapsed_seconds,
      completed, length(tasks), proc.time()[["elapsed"]] - study_started))
    if (!is.null(run_dir)) .gfb_atomic_rds(list(schema_version = settings$schema_version,
      completed = completed, total = length(tasks), completed_indices = which(!vapply(records, is.null, logical(1))),
      partial_dir = paths$partial, statuses = vapply(records, function(r) if (is.null(r)) "pending" else r$status, character(1))),
      file.path(run_dir, "checkpoint.rds"), overwrite = TRUE)
  }
  report_workers <- function(info) {
    for (worker in info) {
      log_line(sprintf("Worker pid=%s; BLAS=%s; fit=%s (%s); inference=%s (%s)",
        worker$pid, worker$BLAS, worker$runtime$fit$backend, worker$runtime$fit$status,
        worker$runtime$inference$backend, worker$runtime$inference$status))
      for (kind in c("fit", "inference")) if (!is.null(worker$runtime[[kind]]$reason))
        log_line(paste("Worker", worker$pid, kind, worker$runtime[[kind]]$reason))
    }
  }
  if (actual_worker_backend == "in_process") {
    # Restore caller thread limits after in-process execution.
    env_names <- c("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")
    old_env <- Sys.getenv(env_names, unset = NA_character_)
    old_blas <- old_omp <- NULL
    if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
      old_blas <- RhpcBLASctl::blas_get_num_procs(); old_omp <- RhpcBLASctl::omp_get_max_threads()
    }
    on.exit({
      for (nm in env_names) if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(old_env[[nm]]), nm))
      if (!is.null(old_blas)) RhpcBLASctl::blas_set_num_threads(old_blas)
      if (!is.null(old_omp)) RhpcBLASctl::omp_set_num_threads(old_omp)
    }, add = TRUE)
    worker_api <- new.env(parent = globalenv())
    worker_info <- list(.gfb_prepare_worker(fns, task_args, runtime_dir,
      worker_fast_fit, worker_fast_inference, debias_method, paths, inference_workers, worker_api))
    report_workers(worker_info)
    if (!is.null(run_dir)) .gfb_atomic_rds(worker_info, file.path(run_dir, "workers.rds"))
    for (index in seq_along(tasks)) record_result(worker_api$.gfb_worker_run(tasks[[index]]), index)
  } else {
    thread_names <- c("VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS")
    env_names <- c(thread_names, "GFB_WORKER_R_HOME", "GFB_WORKER_BLAS_DIR")
    old_env <- Sys.getenv(env_names, unset = NA_character_)
    restore_env <- function() {
      for (nm in env_names) if (is.na(old_env[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(old_env[[nm]]), nm))
    }
    do.call(Sys.setenv, setNames(as.list(rep("1", length(thread_names))), thread_names))
    cluster_args <- list(names = n_workers)
    if (use_accelerate_launcher) {
      Sys.setenv(GFB_WORKER_R_HOME = normalizePath(R.home(), mustWork = TRUE),
                 GFB_WORKER_BLAS_DIR = normalizePath(blas_dir, mustWork = TRUE))
      cluster_args$rscript <- launcher
      log_line("Worker launcher: local macOS Accelerate (one thread per process).")
    } else {
      log_line("Worker launcher: standard Rscript; BLAS details are recorded in workers.rds.")
    }
    cl <- tryCatch(do.call(parallel::makePSOCKcluster, cluster_args), finally = restore_env())
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # Set one-thread BLAS limits before starting workers.
    initialize_worker <- .gfb_prepare_worker
    environment(initialize_worker) <- globalenv()
    worker_info <- parallel::clusterCall(cl, initialize_worker, fns = fns,
      task_args = task_args, runtime_dir = runtime_dir, use_fast_fit = worker_fast_fit,
      use_fast_inference = worker_fast_inference, debias_method = debias_method,
      paths = paths, inference_workers = inference_workers)
    report_workers(worker_info)
    if (!is.null(run_dir)) .gfb_atomic_rds(worker_info, file.path(run_dir, "workers.rds"))
    # Save progress after each completed task.
    send_call <- getFromNamespace("sendCall", "parallel")
    receive <- getFromNamespace("recvOneResult", "parallel")
    run_task <- function(task) .gfb_worker_run(task)
    environment(run_task) <- globalenv()
    for (index in seq_len(n_workers)) send_call(cl[[index]], run_task, list(tasks[[index]]), tag = index)
    last_heartbeat <- proc.time()[["elapsed"]]
    connections <- lapply(cl, function(node) node$con)
    for (done in seq_along(tasks)) {
      while (!any(socketSelect(connections, timeout = 1))) {
        now <- proc.time()[["elapsed"]]
        if (now - last_heartbeat >= 15) {
          progress_files <- list.files(paths$progress, pattern = "^progress_[0-9]+[.]rds$", full.names = TRUE)
          phases <- vapply(progress_files, function(path) {
            state <- tryCatch(readRDS(path), error = function(e) NULL)
            if (is.null(state$phase)) "unknown" else state$phase
          }, character(1))
          counts <- table(phases)
          log_line(sprintf("Progress: completed=%d/%d; wall=%.1fs; %s", completed,
            length(tasks), now - study_started,
            paste(paste0(names(counts), "=", as.integer(counts)), collapse = ", ")))
          last_heartbeat <- now
        }
      }
      received <- receive(cl)
      record_result(received$value, received$tag)
      next_index <- done + n_workers
      if (next_index <= length(tasks)) send_call(cl[[received$node]], run_task,
        list(tasks[[next_index]]), tag = next_index)
    }
  }
  settings$worker_runtime <- worker_info
  settings$elapsed_seconds <- unname(proc.time()[["elapsed"]] - study_started)
  if (!is.null(run_dir)) .gfb_atomic_rds(settings, file.path(run_dir, "settings.rds"), overwrite = TRUE)
  summaries <- .gfb_sim_summaries(records, K)
  result <- c(list(schema_version = settings$schema_version, settings = settings,
    truth = truth, records = records, run_dir = run_dir), summaries)
  class(result) <- c("gfb_simulation", "list")
  if (!is.null(run_dir)) {
    .gfb_atomic_rds(result, file.path(run_dir, "results.rds"))
    for (nm in names(summaries)) utils::write.csv(summaries[[nm]],
      file.path(run_dir, paste0(nm, ".csv")), row.names = FALSE)
  }
  log_line(sprintf("Finished %d/%d tasks in %.2fs. Failed tasks are retained.",
    completed, length(tasks), settings$elapsed_seconds))
  result
}
