# Check configuration helpers without running the simulation notebook.
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
if (length(script)) script <- gsub("~+~", " ", script, fixed = TRUE)
root <- if (length(script)) dirname(dirname(normalizePath(script[[1L]]))) else getwd()
source_file <- file.path(root, "simulations", "node_inference.Rmd")
stopifnot(file.exists(source_file))
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAIL: ", label, call. = FALSE)
  checks <<- checks + 1L
}
equal <- function(x, y, label) {
  result <- all.equal(x, y, check.attributes = FALSE)
  if (!isTRUE(result)) {
    stop("FAIL: ", label, ": ", paste(result, collapse = "; "), call. = FALSE)
  }
  checks <<- checks + 1L
}
error <- function(expr, pattern) {
  caught <- tryCatch({ force(expr); NULL }, error = identity)
  check(inherits(caught, "error") && grepl(pattern, conditionMessage(caught)),
        paste("expected error matching", pattern))
}

lines <- readLines(source_file, warn = FALSE)
in_chunk <- FALSE
r_code <- character()
for (line in lines) {
  if (!in_chunk && grepl("^```\\{r", line)) { in_chunk <- TRUE; next }
  if (in_chunk && grepl("^```\\s*$", line)) { in_chunk <- FALSE; next }
  if (in_chunk) r_code <- c(r_code, line)
}
check(!in_chunk, "R code chunks are closed")
expressions <- parse(text = r_code, keep.source = FALSE)
assignment_to <- function(expr, name = NULL) {
  is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
    is.symbol(expr[[2L]]) &&
    (is.null(name) || identical(as.character(expr[[2L]]), name))
}
needed <- c("normalize_model_name", "parse_lambda_values",
            "resolve_latent_parameters", "default_model_configs",
            "get_model_config", "validate_model_config", "make_param_overrides",
            "number_tag", "numeric_vector_tag", "simulation_output_dir")
api <- new.env(parent = baseenv())
for (expr in expressions) {
  if (assignment_to(expr) && as.character(expr[[2L]]) %in% needed &&
      is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function"))) {
    eval(expr, api)
  }
}
check(setequal(ls(api), needed), "all configuration helpers loaded without runtime code")
check(all(vapply(as.list(api), is.function, logical(1L))),
      "extraction evaluates function definitions only")

# Inspect fallback values without evaluating top-level code.
fallback <- Filter(function(expr) assignment_to(expr, "default_runtime_params"),
                   as.list(expressions))
check(length(fallback) == 1L, "one runtime fallback definition")
fallback_rhs <- fallback[[1L]][[3L]]
check(is.call(fallback_rhs) && identical(fallback_rhs[[1L]], as.name("list")),
      "runtime fallback is a literal list")
runtime <- lapply(as.list(fallback_rhs)[-1L], function(value) {
  if (!is.null(value) && (!is.atomic(value) || length(value) != 1L)) {
    stop("Runtime defaults must be literal scalar values or NULL.", call. = FALSE)
  }
  value
})

if (!requireNamespace("knitr", quietly = TRUE)) {
  stop("Package 'knitr' is required to inspect the notebook header.", call. = FALSE)
}
header_end <- which(lines == "---")
check(length(header_end) >= 2L && header_end[[1L]] == 1L,
      "notebook begins with a YAML header")
# Parse YAML like knitr, preserving the key n and disabling R evaluation.
yaml_params <- lapply(knitr::knit_params(lines, evaluate = FALSE), `[[`, "value")
check(is.list(yaml_params), "YAML supplies simulation parameters")
check(identical(yaml_params$model, "all") && isTRUE(yaml_params$run_simulation),
      "knitting requests all four models and enables simulation")
check(identical(runtime$model, "all") && identical(runtime$run_simulation, FALSE),
      "sourced fallback requests all models without running Monte Carlo")
for (params in list(yaml_params, runtime)) {
  check(is.null(params$k), "K is inferred unless explicitly supplied")
  check(identical(params$truth_check_action, "warning"),
        "effective truth diagnostics warn by default")
  equal(params$n_cores, 7, "default worker cap is seven")
  check("output_dir" %in% names(params) && is.null(params$output_dir),
        "output directory is automatic by default")
}
# Sample size and replication count are user-editable, not fixed test targets.
for (field in c("n", "n_sim")) {
  value <- yaml_params[[field]]
  check(is.null(value) || (is.numeric(value) && length(value) == 1L &&
        is.finite(value) && value >= 1 && value == floor(value)),
        paste("YAML", field, "is a positive integer or NULL"))
}
equal(api$parse_lambda_values(yaml_params$lambda_values), c(0.5, 0.5),
      "YAML Lambda diagonal")
equal(yaml_params$target_node, 1, "YAML target node")

models <- c("gaussian", "poisson", "bernoulli_logit", "bernoulli_probit")
paper <- list(n = 1000L, k = 2L, n_sim = 1000L, lambda_vec = c(0.5, 0.5),
              target_node = 1L, design_seed = 238687L, replication_seed_base = 2000L)
optimizer <- list(
  gaussian = c(max_iter = 10000, learning_rate = 0.1, tolerance = 1e-3),
  poisson = c(max_iter = 60000, learning_rate = 0.02, tolerance = 1e-3),
  bernoulli_logit = c(max_iter = 40000, learning_rate = 1, tolerance = 1e-3),
  bernoulli_probit = c(max_iter = 40000, learning_rate = 1, tolerance = 1e-3)
)
config_from_params <- function(model, params) {
  base <- api$default_model_configs()[[model]]
  api$get_model_config(model, api$make_param_overrides(base, params))
}
for (model in models) {
  config <- api$get_model_config(model)
  equal(config[names(paper)], paper, paste(model, "bare defaults are valid paper settings"))
  equal(unlist(config[names(optimizer[[model]])]), optimizer[[model]],
        paste(model, "retains model-specific optimizer settings"))
  for (mode in c("render", "source")) {
    params <- if (mode == "render") yaml_params else runtime
    effective <- config_from_params(model, params)
    expected <- paper
    if (mode == "render") {
      for (field in c("n", "n_sim")) {
        if (!is.null(yaml_params[[field]])) expected[[field]] <- as.integer(yaml_params[[field]])
      }
    }
    equal(effective[names(paper)], expected, paste(model, mode, "selected experiment settings"))
    check(identical(effective$truth_check_action, "warning"),
          paste(model, "render/source truth warning override"))
  }
}

# Both trial and paper sizes must override the same model defaults correctly.
for (model in models) for (size in c(100L, 1000L)) {
  selected <- yaml_params
  selected$n <- size
  selected$n_sim <- size
  effective <- config_from_params(model, selected)
  equal(effective[c("n", "n_sim")], list(n = size, n_sim = size),
        paste(model, "explicit sample size and replication override", size))
}

# Dimension inference and explicit K must agree for supported experiment sizes.
lambda_cases <- list(0.7, c(0.5, -0.35), c(0.55, 0.4, 0.2))
for (model in models) {
  for (k in seq_along(lambda_cases)) {
    lambda <- lambda_cases[[k]]
    params <- list(n = 24, n_sim = 3, k = NULL,
                   lambda_values = paste(lambda, collapse = ","), target_node = 24)
    inferred <- config_from_params(model, params)
    params$k <- k
    explicit <- config_from_params(model, params)
    check(identical(inferred, explicit), paste(model, "K", k, "explicit/inferred parity"))
    equal(inferred[c("n", "k", "n_sim", "lambda_vec", "target_node")],
          list(n = 24L, k = as.integer(k), n_sim = 3L,
               lambda_vec = lambda, target_node = 24L),
          paste(model, "K", k, "overrides take precedence"))
  }
}
base <- api$get_model_config("gaussian")
for (bad_k in list(0, -1, 1.5, NA_real_, Inf, c(1, 2))) {
  error(api$resolve_latent_parameters(base, bad_k, "0.5,0.5"),
        "k must be one positive integer")
}
error(config_from_params("gaussian", list(k = 3, lambda_values = "0.5,0.5")),
      "k=3 but lambda_values contains 2 entries")
error(config_from_params("gaussian", list(k = 3, lambda_values = NULL)),
      "provide lambda_values")
error(api$get_model_config("gaussian", list(n = 2)), "Require n > k")
for (bad_node in list(0, -1, 1001, 1.5, NA_real_, c(1, 2))) {
  error(api$get_model_config("gaussian", list(target_node = bad_node)), "target_node")
}
for (lambda in list(c(0.5), c(0.5, NA_real_), c(0.5, 0))) {
  error(api$get_model_config("gaussian", list(lambda_vec = lambda)),
        "lambda_vec|lambda value")
}

# Seed overrides stay separate; validate the entire replication seed range.
changed_seeds <- config_from_params("gaussian", list(
  design_seed = 12345, replication_seed_base = 9876))
equal(changed_seeds[c("design_seed", "replication_seed_base")],
      list(design_seed = 12345L, replication_seed_base = 9876L), "seed overrides")
equal(api$get_model_config("gaussian", list(replication_seed_base = 9876))$design_seed,
      paper$design_seed, "replication seed override preserves latent design seed")
for (bad in list(list(design_seed = -1), list(design_seed = 1.5),
                 list(replication_seed_base = -1),
                 list(replication_seed_base = .Machine$integer.max - 999L))) {
  error(api$get_model_config("gaussian", bad), "Seeds|integer scalars")
}
edge_seed <- api$get_model_config("gaussian", list(
  design_seed = .Machine$integer.max,
  replication_seed_base = .Machine$integer.max - 1000L))
equal(edge_seed$replication_seed_base + edge_seed$n_sim, .Machine$integer.max,
      "largest legal final replication seed is accepted")

# Automatic folders encode experiment settings and always separate models.
default_path <- api$simulation_output_dir(base)
expected_parent <- file.path("..", "data", "node1", paste0(
  "n1000_k2_lambda_", api$numeric_vector_tag(c(0.5, 0.5)), "_target1"))
equal(default_path, file.path(expected_parent, "gaussian"), "automatic paper folder")
for (model in models) {
  config <- api$get_model_config(model)
  expected <- file.path(expected_parent, sub("^bernoulli_", "", model))
  equal(api$simulation_output_dir(config), expected, paste(model, "automatic folder"))
  equal(api$simulation_output_dir(config, multiple_models = TRUE), expected,
        paste(model, "automatic folder is identical for single/all runs"))
  equal(api$simulation_output_dir(config, "custom results", multiple_models = FALSE),
        "custom results", paste(model, "explicit single-model path is exact"))
  equal(api$simulation_output_dir(config, "custom results", multiple_models = TRUE),
        file.path("custom results", sub("^bernoulli_", "", model)),
        paste(model, "explicit shared path receives model subfolder"))
}
variations <- list(list(n = 500), list(lambda_vec = c(0.5, 0.4)),
                   list(lambda_vec = c(0.5, 0.5000000001)), list(target_node = 2),
                   list(k = 1, lambda_vec = 0.5))
paths <- vapply(variations, function(overrides) {
  api$simulation_output_dir(api$get_model_config("gaussian", overrides))
}, character(1L))
check(length(unique(c(default_path, paths))) == length(paths) + 1L,
      "n, Lambda (including rounding collisions), target, and K distinguish paths")
cat(sprintf("PASS: %d simulation configuration checks; no truth generation or Monte Carlo was run.\n",
            checks))
