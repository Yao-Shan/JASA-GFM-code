# GFM estimation for one observed network (base R only).
# Public interface: gfm_read_matrix(), gfm_fit(), gfm_predict().
# See practitioner/README.txt for input conventions and interpretation.
# Core formulas extracted from code/node_inference.Rmd on 2026-10-12.
# Original simulation files are unchanged. Sourcing only defines functions.
# The plug-in correction uses U_anchor = U_hat in the fitted coordinate frame.
# Practical additions: validated inputs, explicit failures, resource guard,
# and rank checks before accepting the correction. Defaults/formulas match
# the original estimator on regular inputs. All functions use base R.

# Resolve the four observation-model names and accepted aliases.
.gfm_normalize_model_name <- function(model) {
  key <- tolower(gsub("[- ]", "_", model))
  aliases <- c(
    logit = "bernoulli_logit",
    bernoulli_logit = "bernoulli_logit",
    probit = "bernoulli_probit",
    bernoulli_probit = "bernoulli_probit",
    gaussian = "gaussian",
    normal = "gaussian",
    poisson = "poisson"
  )
  if (!key %in% names(aliases)) {
    stop(
      "Unknown model: ", model,
      ". Choose bernoulli_logit, bernoulli_probit, gaussian, or poisson."
    )
  }
  unname(aliases[[key]])
}

# Remove numerical asymmetry from a square matrix.
.gfm_symmetrize <- function(A) (A + t(A)) / 2

# Evaluate the logistic mean without overflow.
.gfm_sigmoid <- function(x) {
  out <- numeric(length(x))
  positive <- x >= 0
  out[positive] <- 1 / (1 + exp(-x[positive]))
  exp_x <- exp(x[!positive])
  out[!positive] <- exp_x / (1 + exp_x)
  dim(out) <- dim(x)
  out
}

# Evaluate the logistic log normalizer stably.
.gfm_stable_log1pexp <- function(x) {
  pmax(x, 0) + log1p(exp(-abs(x)))
}

# Invert positive information eigenvalues after removing constraint directions.
.gfm_safe_pinv_psd <- function(
  A, n_drop = 0L, rel_tol = 1e-10, negative_tol = 1e-7,
  require_rank = FALSE
) {
  A <- .gfm_symmetrize(A)
  eg <- eigen(A, symmetric = TRUE)
  vals <- eg$values
  p <- length(vals)
  if (any(!is.finite(vals))) stop("Non-finite information eigenvalue.")
  if (n_drop < 0L || n_drop > p) stop("Invalid n_drop.")
  if (p == 0L || n_drop == p) return(matrix(0, p, p))

  scale <- max(abs(vals))
  if (scale == 0) {
    if (require_rank) stop("Projected node information is zero.")
    return(matrix(0, p, p))
  }
  if (min(vals) < -negative_tol * scale) {
    stop(sprintf(
      "Matrix expected to be PSD has min eigenvalue %.4e (scale %.4e).",
      min(vals), scale
    ))
  }

  keep <- rep(TRUE, p)
  if (n_drop > 0L) {
    keep[order(abs(vals))[seq_len(n_drop)]] <- FALSE
  }
  keep <- keep & vals > rel_tol * scale
  if (require_rank && sum(keep) != p - n_drop) {
    stop("Projected node information has extra unidentified directions.")
  }
  vals_inv <- numeric(p)
  vals_inv[keep] <- 1 / vals[keep]
  .gfm_symmetrize(eg$vectors %*% diag(vals_inv, p) %*% t(eg$vectors))
}

# Reject singular or ill-conditioned information matrices.
.gfm_assert_positive_definite <- function(
  A, label = "Information matrix", rel_tol = 1e-10
) {
  values <- eigen(
    .gfm_symmetrize(A), symmetric = TRUE, only.values = TRUE
  )$values
  if (!length(values) || any(!is.finite(values))) {
    stop(label, " has non-finite or missing eigenvalues.")
  }
  scale <- max(abs(values))
  if (scale == 0 || min(values) <= rel_tol * scale) {
    stop(sprintf(
      "%s is not numerically positive definite: min %.4e, max scale %.4e.",
      label, min(values), scale
    ))
  }
  invisible(values)
}

# Invert a positive definite matrix and report numerical diagnostics.
.gfm_ordinary_inverse_pd <- function(
  A, label = "Information matrix", rel_tol = 1e-10
) {
  A <- .gfm_symmetrize(A)
  values <- .gfm_assert_positive_definite(A, label, rel_tol)
  inverse <- .gfm_symmetrize(chol2inv(chol(A)))
  inverse_residual <- max(abs(A %*% inverse - diag(nrow(A))))
  list(
    inverse = inverse,
    eigenvalues = values,
    min_eigenvalue = min(values),
    max_eigenvalue = max(values),
    condition_number = max(values) / min(values),
    inverse_residual = inverse_residual
  )
}

# Vectorize Lambda using diagonal/sqrt(2), then row-wise upper pairs.
# The result has d = K(K+1)/2 entries; this is a parameter vector, not a
# standardized error. For K=2: (Lambda_11/sqrt(2), Lambda_22/sqrt(2), Lambda_12).
.gfm_get_lambda_v <- function(Lambda) {
  k <- nrow(Lambda)
  out <- diag(Lambda) / sqrt(2)
  if (k >= 2L) {
    for (a in seq_len(k - 1L)) {
      for (b in seq.int(a + 1L, k)) out <- c(out, Lambda[a, b])
    }
  }
  as.vector(out)
}

# Differentiate one edge predictor with respect to lambda_v.
.gfm_lambda_direction_derivative <- function(ui, uj) {
  k <- length(ui)
  out <- sqrt(2) * ui * uj
  if (k >= 2L) {
    for (a in seq_len(k - 1L)) {
      for (b in seq.int(a + 1L, k)) {
        out <- c(out, ui[a] * uj[b] + ui[b] * uj[a])
      }
    }
  }
  as.vector(out)
}

# Project onto U' U = n I and certify the row bound; never row-clip.
.gfm_project_C_U_certified <- function(
  Z, bound, row_tolerance = 1e-10, gram_tolerance = 1e-10,
  objective_tolerance = 1e-10, context = ""
) {
  if (!is.matrix(Z) || !is.numeric(Z) || any(!is.finite(Z))) {
    stop("Z must be a finite numeric matrix in the C_U projection.")
  }
  n <- nrow(Z)
  k <- ncol(Z)
  if (k < 1L || n < k || length(bound) != 1L ||
      !is.finite(bound) || bound < sqrt(k)) {
    stop("C_U is infeasible: require n >= k and M_bound_U >= sqrt(k).")
  }

  # The scaled polar factor is a global Euclidean projection onto
  # {U: U^T U = n I_k}. If it also satisfies the row-norm bound, it is
  # therefore a global projection onto the intersection C_U. When the row
  # bound is active, clipping would destroy U^T U = n I_k; fail explicitly
  # instead of silently replacing the theoretical joint projection.
  decomposition <- svd(Z, nu = k, nv = k)
  U <- sqrt(n) * decomposition$u[, seq_len(k), drop = FALSE] %*%
    t(decomposition$v[, seq_len(k), drop = FALSE])
  gram_error <- max(abs(crossprod(U) / n - diag(k)))
  max_row_norm <- max(sqrt(rowSums(U^2)))
  row_violation <- max(0, max_row_norm - bound)
  objective <- 0.5 * sum((U - Z)^2)
  stiefel_lower_bound <- 0.5 * (
    sum(Z^2) + n * k - 2 * sqrt(n) * sum(decomposition$d)
  )
  objective_gap <- abs(objective - stiefel_lower_bound)
  objective_scale <- max(1, abs(objective), abs(stiefel_lower_bound))
  if (gram_error > gram_tolerance) {
    stop(sprintf(
      "The scaled-polar projection failed its Gram check: %.4e > %.4e.",
      gram_error, gram_tolerance
    ))
  }
  if (objective_gap > objective_tolerance * objective_scale) {
    stop(sprintf(
      "The scaled-polar objective certificate failed: gap %.4e.",
      objective_gap
    ))
  }
  if (row_violation > row_tolerance * max(1, bound)) {
    message_text <- sprintf(
      paste0(
        "Certified joint C_U projection unavailable%s: scaled-polar row ",
        "norm %.16g exceeds M_bound_U %.16g by %.4e. The run must stop; ",
        "row clipping is not a joint Euclidean projection."
      ),
      if (nzchar(context)) paste0(" [", context, "]") else "",
      max_row_norm, bound, row_violation
    )
    condition <- structure(
      list(
        message = message_text, call = NULL, context = context,
        row_max = max_row_norm, bound = bound,
        row_violation = row_violation, gram_error = gram_error,
        objective_gap = objective_gap
      ),
      class = c("active_U_row_bound", "error", "condition")
    )
    stop(condition)
  }

  list(
    U = U,
    method = "polar_exact",
    max_row_norm = max_row_norm,
    max_row_ratio = max_row_norm / bound,
    row_violation = row_violation,
    row_margin = bound - max_row_norm,
    gram_error = gram_error,
    objective_gap = objective_gap
  )
}

# Project symmetric Lambda onto the spectral-norm bound.
.gfm_proximal_projection_Lambda <- function(Lambda, bound) {
  eg <- eigen(.gfm_symmetrize(Lambda), symmetric = TRUE)
  vals <- pmax(pmin(eg$values, bound), -bound)
  .gfm_symmetrize(eg$vectors %*% diag(vals, length(vals)) %*% t(eg$vectors))
}

# Measure constrained stationarity with the original effective step sizes.
.gfm_projected_gradient_mapping_components <- function(
  U, U_next, Lambda, Lambda_next, learning_rate, n
) {
  if (!is.finite(learning_rate) || learning_rate <= 0 ||
      !is.finite(n) || n <= 0) {
    stop("learning_rate and n must be positive in the residual calculation.")
  }
  eta_U <- learning_rate / n
  eta_Lambda_v <- learning_rate / (2 * n^2)
  c(
    U = max(abs(U_next - U)) / eta_U,
    Lambda_v = max(abs(
      .gfm_get_lambda_v(Lambda_next) - .gfm_get_lambda_v(Lambda)
    )) / eta_Lambda_v
  )
}

# Compute stable probit Mills ratios using log probabilities.
.gfm_probit_ratios <- function(theta) {
  log_phi <- dnorm(theta, log = TRUE)
  log_p1 <- pnorm(theta, log.p = TRUE)
  log_p0 <- pnorm(theta, lower.tail = FALSE, log.p = TRUE)
  list(
    p1 = pnorm(theta),
    p0 = pnorm(theta, lower.tail = FALSE),
    a = exp(log_phi - log_p1),
    b = exp(log_phi - log_p0),
    log_phi = log_phi,
    log_p1 = log_p1,
    log_p0 = log_p0
  )
}

# Map linear predictors to edge means.
.gfm_link_probability <- function(theta, config) {
  switch(
    config$model,
    bernoulli_logit = .gfm_sigmoid(theta),
    bernoulli_probit = pnorm(theta),
    poisson = exp(theta),
    gaussian = theta
  )
}

# Map the spectral initializer to the model's predictor scale.
.gfm_link_inverse <- function(value, config, eps = config$probability_eps) {
  switch(
    config$model,
    bernoulli_logit = {
      p <- pmax(pmin(value, 1 - eps), eps)
      log(p / (1 - p))
    },
    bernoulli_probit = qnorm(pmax(pmin(value, 1 - eps), eps)),
    poisson = log(pmax(value, 0.1)),
    gaussian = value
  )
}

# Evaluate the off-diagonal negative log likelihood up to constants.
.gfm_negative_loglikelihood_edges <- function(theta, y, config) {
  switch(
    config$model,
    bernoulli_logit = sum(.gfm_stable_log1pexp(theta) - y * theta),
    bernoulli_probit = {
      log_p1 <- pnorm(theta, log.p = TRUE)
      log_p0 <- pnorm(theta, lower.tail = FALSE, log.p = TRUE)
      -sum(ifelse(y > 0.5, log_p1, log_p0))
    },
    gaussian = sum((theta - y)^2) / (2 * config$sigma^2),
    poisson = sum(exp(theta) - y * theta)
  )
}

# Evaluate edgewise negative likelihood derivatives.
.gfm_negative_score_matrix <- function(theta, A, config) {
  switch(
    config$model,
    bernoulli_logit = .gfm_sigmoid(theta) - A,
    bernoulli_probit = {
      r <- .gfm_probit_ratios(theta)
      out <- ifelse(A > 0.5, -r$a, r$b)
      dim(out) <- dim(theta)
      out
    },
    gaussian = (theta - A) / config$sigma^2,
    poisson = exp(theta) - A
  )
}

# Compute Fisher information and expected higher derivatives for bias terms.
# Here ell is one edge's log likelihood: information = E[-ell''],
# expected_ell3 = E[ell'''], and expected_score_ell2 = E[ell' * ell''].
.gfm_expected_derivative_terms <- function(theta, config) {
  zero <- 0 * theta

  if (config$model == "bernoulli_logit") {
    p <- .gfm_sigmoid(theta)
    information <- p * (1 - p)
    expected_ell3 <- -information * (1 - 2 * p)
    return(list(
      mean_parameter = p,
      information = information,
      expected_ell3 = expected_ell3,
      expected_score_ell2 = zero,
      observed_information_y1 = information,
      observed_information_y0 = information
    ))
  }

  if (config$model == "bernoulli_probit") {
    r <- .gfm_probit_ratios(theta)
    a <- r$a
    b <- r$b
    ell2_1 <- -theta * a - a^2
    ell3_1 <- a * (theta^2 - 1) + 3 * theta * a^2 + 2 * a^3
    ell2_0 <- theta * b - b^2
    ell3_0 <- b * (1 - theta^2) + 3 * theta * b^2 - 2 * b^3
    information <- exp(2 * r$log_phi - r$log_p1 - r$log_p0)
    return(list(
      mean_parameter = r$p1,
      information = information,
      expected_ell3 = information * (3 * theta + 2 * a - 2 * b),
      expected_score_ell2 = information * (b - a - theta),
      observed_information_y1 = -ell2_1,
      observed_information_y0 = -ell2_0,
      ell3_y1 = ell3_1,
      ell3_y0 = ell3_0
    ))
  }

  if (config$model == "gaussian") {
    return(list(
      mean_parameter = theta,
      information = zero + 1 / config$sigma^2,
      expected_ell3 = zero,
      expected_score_ell2 = zero,
      observed_information_y1 = zero + 1 / config$sigma^2,
      observed_information_y0 = zero + 1 / config$sigma^2
    ))
  }

  intensity <- exp(theta)
  list(
    mean_parameter = intensity,
    information = intensity,
    expected_ell3 = -intensity,
    expected_score_ell2 = zero,
    observed_information_y1 = intensity,
    observed_information_y0 = intensity
  )
}

# Compute realized edge information, including outcome dependence for probit.
.gfm_observed_information_matrix <- function(theta, config, A = NULL) {
  # For Gaussian, logit, and Poisson, the observed negative Hessian is
  # identical to the expected Fisher information.  Probit is different:
  # its observed curvature depends on the realized binary outcome.
  if (config$model != "bernoulli_probit") {
    return(.gfm_expected_derivative_terms(theta, config)$information)
  }

  if (is.null(A)) {
    stop("A is required for the observed probit information matrix.")
  }
  if (!identical(dim(A), dim(theta))) {
    stop("A and theta must have the same dimensions.")
  }

  r <- .gfm_probit_ratios(theta)
  information_y1 <- r$a * (theta + r$a)
  information_y0 <- r$b * (r$b - theta)
  information <- ifelse(A > 0.5, information_y1, information_y0)
  dim(information) <- dim(theta)

  # Both expressions are nonnegative analytically.  Remove only negative
  # roundoff caused by tail-ratio cancellation.
  pmax(information, 0)
}

# Select the spectral initializer's retained singular values.
.gfm_spectral_keep_indices <- function(singular_values, k, offset) {
  p <- length(singular_values)
  threshold_index <- k + offset + 1L
  threshold <- if (threshold_index <= p) {
    singular_values[threshold_index] + 1e-6
  } else {
    0
  }
  keep <- which(singular_values > threshold)
  if (length(keep) < k) keep <- seq_len(min(k, p))
  keep
}

# Initialize U and Lambda from the observed network.
.gfm_initialize_spectral <- function(A, config) {
  n <- nrow(A)
  k <- config$k
  sv <- svd(A)

  if (config$model %in% c("bernoulli_logit", "bernoulli_probit")) {
    keep <- .gfm_spectral_keep_indices(sv$d, k, offset = 2L)
    P_tilde <- sv$u[, keep, drop = FALSE] %*%
      diag(sv$d[keep], length(keep)) %*% t(sv$v[, keep, drop = FALSE])
    P_tilde <- pmax(
      pmin(.gfm_symmetrize(P_tilde), 1 - config$probability_eps),
      config$probability_eps
    )
    Theta_tilde <- .gfm_symmetrize(.gfm_link_inverse(P_tilde, config))
  } else if (config$model == "gaussian") {
    keep <- .gfm_spectral_keep_indices(sv$d, k, offset = 0L)
    Theta_tilde <- sv$u[, keep, drop = FALSE] %*%
      diag(sv$d[keep], length(keep)) %*% t(sv$v[, keep, drop = FALSE])
    Theta_tilde <- .gfm_symmetrize(Theta_tilde)
  } else {
    keep <- .gfm_spectral_keep_indices(sv$d, k, offset = 2L)
    intensity_tilde <- sv$u[, keep, drop = FALSE] %*%
      diag(sv$d[keep], length(keep)) %*% t(sv$v[, keep, drop = FALSE])
    intensity_tilde <- pmax(.gfm_symmetrize(intensity_tilde), 0.1)
    Theta_tilde <- .gfm_symmetrize(log(intensity_tilde))
  }

  diag(Theta_tilde) <- 0
  eg <- eigen(Theta_tilde, symmetric = TRUE)
  idx <- order(abs(eg$values), decreasing = TRUE)[seq_len(k)]
  list(
    U_init = sqrt(n) * eg$vectors[, idx, drop = FALSE],
    Lambda_init = diag(eg$values[idx] / n, k, k)
  )
}

# Evaluate the objective on the unordered observed edges.
.gfm_calculate_objective_network <- function(U, Lambda, A, config) {
  Theta <- U %*% Lambda %*% t(U)
  diag(Theta) <- 0
  upper <- upper.tri(A, diag = FALSE)
  .gfm_negative_loglikelihood_edges(Theta[upper], A[upper], config)
}

# Compute the joint U and symmetric-Lambda gradients.
.gfm_calculate_gradient_network <- function(U, Lambda, A, config) {
  Theta <- U %*% Lambda %*% t(U)
  diag(Theta) <- 0
  D <- .gfm_negative_score_matrix(Theta, A, config)
  diag(D) <- 0
  grad_Lambda <- 0.5 * t(U) %*% D %*% U
  list(
    grad_U = D %*% U %*% Lambda,
    grad_Lambda = .gfm_symmetrize(grad_Lambda)
  )
}

# Fit one network with the original fixed-step projected gradient algorithm.
# Return raw matrices and diagnostics. Codes: 0 = stopping criterion met,
# 1 = iteration limit reached, 2 = non-finite residual or objective.
.gfm_run_pgd_network <- function(A, config) {
  n <- config$n
  init <- .gfm_initialize_spectral(A, config)
  U_projection <- .gfm_project_C_U_certified(
    init$U_init, config$M_bound_U, context = "spectral initialization"
  )
  U_hat <- U_projection$U
  Lambda_hat <- .gfm_proximal_projection_Lambda(
    init$Lambda_init, config$M_bound_Lambda
  )

  convergence_code <- 1L
  scaled_gradient <- Inf
  projected_residual_U <- Inf
  projected_residual_Lambda_v <- Inf
  iter_done <- 0L
  projection_calls <- 1L
  projection_max_row_ratio <- U_projection$max_row_ratio
  projection_min_row_margin <- U_projection$row_margin
  projection_max_gram_error <- U_projection$gram_error
  projection_max_objective_gap <- U_projection$objective_gap

  for (iter in seq_len(config$max_iter)) {
    grads <- .gfm_calculate_gradient_network(U_hat, Lambda_hat, A, config)
    U_projection <- .gfm_project_C_U_certified(
      U_hat - config$learning_rate * grads$grad_U / n,
      config$M_bound_U, context = sprintf("PGD iteration %d", iter)
    )
    U_next <- U_projection$U
    projection_calls <- projection_calls + 1L
    projection_max_row_ratio <- max(
      projection_max_row_ratio, U_projection$max_row_ratio
    )
    projection_min_row_margin <- min(
      projection_min_row_margin, U_projection$row_margin
    )
    projection_max_gram_error <- max(
      projection_max_gram_error, U_projection$gram_error
    )
    projection_max_objective_gap <- max(
      projection_max_objective_gap, U_projection$objective_gap
    )
    Lambda_next <- .gfm_proximal_projection_Lambda(
      Lambda_hat - config$learning_rate * grads$grad_Lambda / n^2,
      config$M_bound_Lambda
    )

    # On the orthogonality-constrained parameter space, the raw Euclidean
    # gradient need not vanish at the optimum. The mapping must be divided by
    # the effective steps actually used by the code: eta_U = lr/n and, under
    # the paper's lambda_v parameterization, eta_lambda_v = lr/(2 n^2).
    residuals <- .gfm_projected_gradient_mapping_components(
      U_hat, U_next, Lambda_hat, Lambda_next,
      config$learning_rate, n
    )
    projected_residual_U <- residuals[["U"]]
    projected_residual_Lambda_v <- residuals[["Lambda_v"]]
    scaled_gradient <- max(residuals)
    if (!is.finite(scaled_gradient)) {
      convergence_code <- 2L
      break
    }
    U_hat <- U_next
    Lambda_hat <- Lambda_next
    iter_done <- iter
    if (scaled_gradient < config$tolerance) {
      convergence_code <- 0L
      break
    }
  }

  # The final update can itself cross the tolerance. Re-evaluate at the
  # returned iterate so a max-iteration exit is classified correctly.
  if (convergence_code == 1L) {
    final_grads <- .gfm_calculate_gradient_network(U_hat, Lambda_hat, A, config)
    final_U_projection <- .gfm_project_C_U_certified(
      U_hat - config$learning_rate * final_grads$grad_U / n,
      config$M_bound_U, context = "final residual check"
    )
    final_U_next <- final_U_projection$U
    projection_calls <- projection_calls + 1L
    projection_max_row_ratio <- max(
      projection_max_row_ratio, final_U_projection$max_row_ratio
    )
    projection_min_row_margin <- min(
      projection_min_row_margin, final_U_projection$row_margin
    )
    projection_max_gram_error <- max(
      projection_max_gram_error, final_U_projection$gram_error
    )
    projection_max_objective_gap <- max(
      projection_max_objective_gap, final_U_projection$objective_gap
    )
    final_Lambda_next <- .gfm_proximal_projection_Lambda(
      Lambda_hat -
        config$learning_rate * final_grads$grad_Lambda / n^2,
      config$M_bound_Lambda
    )
    final_residuals <- .gfm_projected_gradient_mapping_components(
      U_hat, final_U_next, Lambda_hat, final_Lambda_next,
      config$learning_rate, n
    )
    projected_residual_U <- final_residuals[["U"]]
    projected_residual_Lambda_v <- final_residuals[["Lambda_v"]]
    scaled_gradient <- max(final_residuals)
    if (!is.finite(scaled_gradient)) {
      convergence_code <- 2L
    } else if (scaled_gradient < config$tolerance) {
      convergence_code <- 0L
    }
  }

  final_objective <- .gfm_calculate_objective_network(U_hat, Lambda_hat, A, config)
  if (!is.finite(final_objective)) convergence_code <- 2L
  list(
    U_hat = U_hat,
    Lambda_hat = .gfm_symmetrize(Lambda_hat),
    convergence_code = convergence_code,
    n_iter = iter_done,
    scaled_gradient = scaled_gradient,
    projected_residual_U = projected_residual_U,
    projected_residual_Lambda_v = projected_residual_Lambda_v,
    projection_calls = projection_calls,
    projection_max_row_ratio = projection_max_row_ratio,
    projection_min_row_margin = projection_min_row_margin,
    projection_max_gram_error = projection_max_gram_error,
    projection_max_objective_gap = projection_max_objective_gap,
    final_objective = final_objective
  )
}

# Build normalization and rotation-constraint derivatives in the chosen frame.
# DUG is (nK)-by-K^2. Rows stack U node by node, not R's column-major order.
# Columns cover K diagonal normalization constraints, K(K-1)/2 cross-products,
# and K(K-1)/2 rotation constraints anchored at U_anchor.
.gfm_build_D_UG <- function(U, U_anchor) {
  n <- nrow(U)
  k <- ncol(U)
  if (!identical(dim(U_anchor), dim(U))) {
    stop("U and U_anchor must have identical dimensions in .gfm_build_D_UG().")
  }
  DUG <- matrix(0, n * k, k^2)

  for (a in seq_len(k)) {
    DUG[seq(a, n * k, by = k), a] <- U[, a]
  }

  pair <- 1L
  if (k >= 2L) {
    for (a in seq_len(k - 1L)) {
      for (b in seq.int(a + 1L, k)) {
        col_g2 <- k + pair
        col_g3 <- k + k * (k - 1L) / 2L + pair
        rows_a <- seq(a, n * k, by = k)
        rows_b <- seq(b, n * k, by = k)
        DUG[rows_a, col_g2] <- U[, b] / sqrt(2)
        DUG[rows_b, col_g2] <- U[, a] / sqrt(2)
        DUG[rows_a, col_g3] <- U_anchor[, b] / sqrt(2)
        DUG[rows_b, col_g3] <- -U_anchor[, a] / sqrt(2)
        pair <- pair + 1L
      }
    }
  }

  if (any(!is.finite(DUG))) stop("D_UG contains non-finite entries.")
  DUG
}

# Assemble observed information and the exact nuisance projector.
# SUU, SUL, and SLL have sizes (nK)x(nK), (nK)xd, and dxd, where
# d = K(K+1)/2. PUU projects onto allowed U directions. Seff is the dxd
# Schur complement after accounting for those directions. Expected KxK
# node blocks are stored separately for the bias formulas.
.gfm_build_information_blocks <- function(U, Lambda, U_anchor, config, A = NULL) {
  n <- nrow(U)
  k <- ncol(U)
  d <- k * (k + 1L) / 2L
  Theta <- U %*% Lambda %*% t(U)
  diag(Theta) <- 0
  # Match the PDF's realized negative-Hessian information blocks. For probit
  # this uses the observed outcomes; for the other models it equals Fisher.
  V <- .gfm_observed_information_matrix(Theta, config, A)
  V_expected <- .gfm_expected_derivative_terms(Theta, config)$information
  diag(V) <- 0
  diag(V_expected) <- 0
  W <- U %*% Lambda

  SUU <- matrix(0, n * k, n * k)
  SUL <- matrix(0, n * k, d)
  SLL <- matrix(0, d, d)
  expected_blocks <- lapply(seq_len(n), function(i) matrix(0, k, k))

  for (i in seq_len(n - 1L)) {
    ii <- ((i - 1L) * k + 1L):(i * k)
    for (j in seq.int(i + 1L, n)) {
      jj <- ((j - 1L) * k + 1L):(j * k)
      v <- V[i, j]
      v_expected <- V_expected[i, j]
      wi <- W[i, ]
      wj <- W[j, ]
      ld <- .gfm_lambda_direction_derivative(U[i, ], U[j, ])

      SUU[ii, ii] <- SUU[ii, ii] + v * tcrossprod(wj)
      SUU[jj, jj] <- SUU[jj, jj] + v * tcrossprod(wi)
      cross <- v * tcrossprod(wj, wi)
      SUU[ii, jj] <- SUU[ii, jj] + cross
      SUU[jj, ii] <- SUU[jj, ii] + t(cross)
      SUL[ii, ] <- SUL[ii, ] + v * tcrossprod(wj, ld)
      SUL[jj, ] <- SUL[jj, ] + v * tcrossprod(wi, ld)
      SLL <- SLL + v * tcrossprod(ld)
      expected_blocks[[i]] <- expected_blocks[[i]] +
        v_expected * tcrossprod(wj)
      expected_blocks[[j]] <- expected_blocks[[j]] +
        v_expected * tcrossprod(wi)
    }
  }

  SUU <- .gfm_symmetrize(SUU)
  SLL <- .gfm_symmetrize(SLL)
  SLU <- t(SUL)
  DUG <- .gfm_build_D_UG(U, U_anchor)
  DGU <- t(DUG)
  DGU_DUG <- .gfm_symmetrize(DGU %*% DUG)
  DGU_DUG_inverse <- .gfm_ordinary_inverse_pd(
    DGU_DUG, "D_GU D_UG nuisance Gram matrix"
  )
  # Finite-sample nuisance projector. Do not replace D_GU D_UG by n I:
  # Pi_nui = D_UG (D_GU D_UG)^+ D_GU. The numerical-rank check above makes
  # the Moore-Penrose inverse equal to the ordinary inverse used here.
  Pi_nui <- .gfm_symmetrize(
    DUG %*% DGU_DUG_inverse$inverse %*% DGU
  )
  PUU <- .gfm_symmetrize(diag(n * k) - Pi_nui)
  PUU_idempotence_error <- max(abs(PUU %*% PUU - PUU))
  PUU_DUG_annihilation_error <- max(abs(PUU %*% DUG))
  if (PUU_idempotence_error > 1e-8 ||
      PUU_DUG_annihilation_error > 1e-8 * max(1, max(abs(DUG)))) {
    stop(sprintf(
      paste0(
        "Failed to construct the finite-sample nuisance projector: ",
        "idempotence error %.4e, annihilation error %.4e."
      ),
      PUU_idempotence_error, PUU_DUG_annihilation_error
    ))
  }
  projected <- .gfm_symmetrize(PUU %*% SUU %*% PUU)
  projected_inv <- .gfm_safe_pinv_psd(
    projected, n_drop = k^2, require_rank = TRUE
  )
  Seff <- .gfm_symmetrize(
    SLL - SLU %*% PUU %*% projected_inv %*% PUU %*% SUL
  )

  blocks <- lapply(seq_len(n), function(i) {
    idx <- ((i - 1L) * k + 1L):(i * k)
    SUU[idx, idx, drop = FALSE]
  })

  list(
    SUU = SUU, SUL = SUL, SLU = SLU, SLL = SLL,
    DUG = DUG, DGU = DGU,
    DGU_DUG = DGU_DUG,
    DGU_DUG_inverse = DGU_DUG_inverse$inverse,
    DUG_gram_error_from_nI = max(abs(DGU_DUG / n - diag(k^2))),
    DUG_gram_condition_number = DGU_DUG_inverse$condition_number,
    Pi_nui = Pi_nui, PUU = PUU,
    PUU_idempotence_error = PUU_idempotence_error,
    PUU_DUG_annihilation_error = PUU_DUG_annihilation_error,
    projected_inv = projected_inv,
    Seff = Seff, blocks = blocks, expected_blocks = expected_blocks
  )
}

# Evaluate the original B1, BM and BP bias formulas with expected node information.
# B1 uses expected third derivatives; BM uses E[ell' * ell'']; BP supplies
# the normalization-related term. U vectors stack nodes and have length nK;
# Lambda vectors use lambda_v and have length d = K(K+1)/2.
.gfm_calculate_bias_terms <- function(U, Lambda, expected_blocks, config) {
  n <- nrow(U)
  k <- ncol(U)
  d <- k * (k + 1L) / 2L
  Theta <- U %*% Lambda %*% t(U)
  diag(Theta) <- 0
  terms <- .gfm_expected_derivative_terms(Theta, config)
  W <- U %*% Lambda
  # A feasible data analysis must report nonidentifiability, not hide it in
  # a zero-eigenvalue pseudoinverse. For regular blocks this is unchanged.
  for (i in seq_along(expected_blocks)) {
    .gfm_assert_positive_definite(
      expected_blocks[[i]], sprintf("Expected information for node %d", i)
    )
  }
  inv_blocks <- lapply(expected_blocks, .gfm_safe_pinv_psd)
  inverse_sum <- Reduce(`+`, inv_blocks)

  B1U <- numeric(n * k)
  B1L <- numeric(d)
  BMU <- numeric(n * k)
  BML <- numeric(d)
  BPU_matrix <- 0.5 * U %*% t(inverse_sum)

  for (i in seq_len(n)) {
    idx_i <- ((i - 1L) * k + 1L):(i * k)
    wi <- matrix(W[i, ], k, 1L)
    b1_i <- matrix(0, k, 1L)
    bm_i <- matrix(0, k, 1L)

    for (j in setdiff(seq_len(n), i)) {
      wj <- matrix(W[j, ], k, 1L)
      ld <- .gfm_lambda_direction_derivative(U[i, ], U[j, ])
      q_i_j <- as.numeric(t(wi) %*% inv_blocks[[j]] %*% wi)
      q_j_i <- as.numeric(t(wj) %*% inv_blocks[[i]] %*% wj)
      e3 <- terms$expected_ell3[i, j]
      em <- terms$expected_score_ell2[i, j]

      b1_i <- b1_i - 0.5 * e3 * (q_i_j + q_j_i) * wj
      bm_i <- bm_i + em * (q_i_j + q_j_i) * wj
      B1L <- B1L - 0.5 * e3 * q_i_j * ld
      BML <- BML + em * q_i_j * ld
    }
    B1U[idx_i] <- as.vector(b1_i)
    BMU[idx_i] <- as.vector(bm_i)
  }

  list(
    B1U = B1U, B1L = B1L,
    BMU = BMU, BML = BML,
    BPU = as.vector(t(BPU_matrix))
  )
}

# Compute additive Lambda corrections with an explicitly supplied coordinate anchor.
# Both returned corrections use lambda_v units. Their sum is ADDED to the
# raw vector; converting to a KxK matrix restores the diagonal sqrt(2) factor.
# U is held fixed here. This function does not construct confidence intervals.
.gfm_calculate_debiasing_components <- function(
  U, Lambda, U_anchor, config, A = NULL
) {
  n <- nrow(U)
  info <- .gfm_build_information_blocks(U, Lambda, U_anchor, config, A = A)
  SLL_inverse <- .gfm_ordinary_inverse_pd(
    info$SLL, "Theorem 2 Lambda information"
  )
  Seff_inverse <- .gfm_ordinary_inverse_pd(
    info$Seff, "Theorem 2 efficient information"
  )
  # The B1/BM/BP formulas use the expected node-information blocks; keep
  # those separate from the observed Hessian blocks used in Seff above.
  bias <- .gfm_calculate_bias_terms(U, Lambda, info$expected_blocks, config)

  core <- (bias$B1L - bias$BML) -
    info$SLU %*% info$PUU %*% info$projected_inv %*% info$PUU %*%
    (bias$B1U - bias$BMU)

  correction_likelihood <- as.vector(
    Seff_inverse$inverse %*% core
  )
  penalty_direction <- info$DUG %*% crossprod(info$DUG, bias$BPU)
  correction_penalty <- -as.vector(
    SLL_inverse$inverse %*% info$SLU %*%
      penalty_direction / n^2
  )

  list(
    Seff = info$Seff,
    lambda_v_hat = .gfm_get_lambda_v(Lambda),
    correction_likelihood = correction_likelihood,
    correction_penalty = correction_penalty,
    correction_total = correction_likelihood + correction_penalty,
    max_abs_BM = max(abs(c(bias$BMU, bias$BML))),
    DUG_gram_error_from_nI = info$DUG_gram_error_from_nI,
    DUG_gram_condition_number = info$DUG_gram_condition_number,
    PUU_idempotence_error = info$PUU_idempotence_error,
    PUU_DUG_annihilation_error = info$PUU_DUG_annihilation_error,
    SLL_min_eigenvalue = SLL_inverse$min_eigenvalue,
    SLL_condition_number = SLL_inverse$condition_number,
    SLL_inverse_residual = SLL_inverse$inverse_residual,
    Seff_min_eigenvalue = Seff_inverse$min_eigenvalue,
    Seff_condition_number = Seff_inverse$condition_number,
    Seff_inverse_residual = Seff_inverse$inverse_residual
  )
}

# ---------------------- Practical public interface ----------------------

# Check one finite scalar; used to reject malformed controls before fitting.
.gfm_positive_scalar <- function(x, name, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 ||
      (integer && (x != floor(x) || x > .Machine$integer.max))) {
    stop(name, " must be a finite positive ",
         if (integer) "integer." else "number.", call. = FALSE)
  }
  invisible(x)
}

# Construct estimation-only controls: no truth, seed, target node or MC count.
.gfm_config <- function(n, K, model, sigma, control) {
  if (!is.character(model) || length(model) != 1L || is.na(model)) {
    stop("model must be one model name.", call. = FALSE)
  }
  model <- .gfm_normalize_model_name(model)
  .gfm_positive_scalar(K, "K", integer = TRUE)
  if (K >= n) stop("K must be smaller than the number of nodes.", call. = FALSE)
  if (model == "gaussian") {
    .gfm_positive_scalar(sigma, "Known Gaussian sigma")
  } else if (!is.null(sigma)) {
    stop("sigma is only used by the Gaussian model; omit it otherwise.", call. = FALSE)
  }
  defaults <- switch(model,
    bernoulli_logit = list(max_iter = 40000L, learning_rate = 1, tolerance = 0.1),
    bernoulli_probit = list(max_iter = 40000L, learning_rate = 1, tolerance = 0.001),
    gaussian = list(max_iter = 10000L, learning_rate = 0.1, tolerance = 0.1),
    poisson = list(max_iter = 60000L, learning_rate = 0.02, tolerance = 0.0001)
  )
  defaults <- c(defaults, list(
    M_bound_U = 20, M_bound_Lambda = 20, probability_eps = 1e-3,
    max_debias_memory_mb = 2048
  ))
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyNA(names(control)) ||
       any(!nzchar(names(control))) || anyDuplicated(names(control))))) {
    stop("control must be a uniquely named list.", call. = FALSE)
  }
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown)) stop("Unknown control: ", paste(unknown, collapse = ", "), call. = FALSE)
  defaults[names(control)] <- control
  for (nm in names(defaults)) {
    .gfm_positive_scalar(defaults[[nm]], paste0("control$", nm), nm == "max_iter")
  }
  if (defaults$M_bound_U < sqrt(K)) {
    stop("control$M_bound_U must be at least sqrt(K).", call. = FALSE)
  }
  if (defaults$probability_eps >= 0.5) {
    stop("control$probability_eps must be below 0.5.", call. = FALSE)
  }
  c(list(n = n, k = as.integer(K), model = model,
         sigma = if (is.null(sigma)) 1 else sigma), defaults)
}

# Validate an undirected complete network; disregard diagonal before all checks.
.gfm_validate_matrix <- function(A, model = NULL) {
  if (!is.matrix(A) || !is.numeric(A) || nrow(A) != ncol(A) || nrow(A) < 2L) {
    stop("A must be a numeric square matrix with at least two nodes.", call. = FALSE)
  }
  rn <- rownames(A)
  cn <- colnames(A)
  if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
    stop("Row and column node names must agree in the same order.", call. = FALSE)
  }
  labels <- if (!is.null(rn)) rn else cn
  if (!is.null(labels) && (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))) {
    stop("Node names must be nonempty and unique.", call. = FALSE)
  }
  diag(A) <- 0
  if (any(!is.finite(A))) stop("All off-diagonal edges must be finite and observed.", call. = FALSE)
  if (!all(A == t(A))) stop("A must be symmetric; edges will not be silently averaged.", call. = FALSE)
  if (!is.null(model)) {
    values <- A[upper.tri(A)]
    if (model %in% c("bernoulli_logit", "bernoulli_probit") &&
        any(values != 0 & values != 1)) {
      stop("Logit and probit require binary (0/1) edges.", call. = FALSE)
    }
    if (model == "poisson" && any(values < 0 | values != floor(values))) {
      stop("Poisson requires nonnegative integer edges.", call. = FALSE)
    }
  }
  if (!is.null(labels)) dimnames(A) <- list(labels, labels)
  A
}

#' Read one observed network from a CSV file.
#'
#' file: CSV path accepted by read.csv(). By default, all cells are data.
#' node_names: TRUE when the first row and first column contain node labels.
#' Returns: a numeric symmetric n-by-n matrix with its diagonal set to zero.
#' Off-diagonal entries must be finite; labels must agree in the same order.
#' This checks the matrix structure; gfm_fit() checks the chosen model's domain.
#' Example: A <- gfm_read_matrix("network.csv")
#' Example with labels: A <- gfm_read_matrix("network.csv", node_names = TRUE)
gfm_read_matrix <- function(file, node_names = FALSE) {
  if (!is.logical(node_names) || length(node_names) != 1L || is.na(node_names)) {
    stop("node_names must be TRUE or FALSE.", call. = FALSE)
  }
  dat <- utils::read.csv(file, header = node_names,
                         row.names = if (node_names) 1L else NULL,
                         check.names = FALSE)
  A <- as.matrix(dat)
  if (!node_names) dimnames(A) <- NULL
  .gfm_validate_matrix(A)
}

# Undo the paper's lambda_v scaling and ordering for arbitrary K.
.gfm_lambda_matrix <- function(v, K) {
  if (length(v) != K * (K + 1L) / 2L) stop("Invalid lambda_v length.")
  Lambda <- diag(sqrt(2) * v[seq_len(K)], K, K)
  pos <- K
  if (K > 1L) {
    for (a in seq_len(K - 1L)) for (b in seq.int(a + 1L, K)) {
      pos <- pos + 1L
      Lambda[a, b] <- Lambda[b, a] <- v[pos]
    }
  }
  Lambda
}

# Label actual parameter entries (not standardized Monte Carlo errors).
.gfm_lambda_names <- function(K) {
  out <- paste0("Lambda[", seq_len(K), ",", seq_len(K), "]/sqrt(2)")
  if (K > 1L) {
    for (a in seq_len(K - 1L)) for (b in seq.int(a + 1L, K)) {
      out <- c(out, paste0("Lambda[", a, ",", b, "]"))
    }
  }
  out
}

#' Estimate a symmetric GFM from one observed network.
#'
#' A: numeric symmetric n-by-n matrix, with all off-diagonal edges observed.
#'    The diagonal is ignored, including NA. Node order is retained.
#' K: positive integer smaller than n; all dimensions adapt to K.
#' model: "logit", "probit", "gaussian", or "poisson" (aliases also accepted).
#'    Bernoulli models require 0/1 edges; Poisson requires nonnegative counts.
#' sigma: known positive Gaussian noise standard deviation; required for
#'    Gaussian and otherwise omitted. The function does not estimate sigma.
#' control: uniquely named list overriding optimizer or resource defaults.
#'    All values are positive finite scalars; max_iter is an integer,
#'    M_bound_U >= sqrt(K), and 0 < probability_eps < 0.5. See README.txt.
#' debias: whether to attempt the additive Lambda correction after convergence.
#'
#' Returns: a gfm_fit list containing U_hat (n-by-K), Lambda_hat (K-by-K),
#'    and, on correction success, Lambda_debiased (K-by-K). These are actual
#'    parameter estimates, not standardized Monte Carlo errors. lambda_v
#'    outputs have K(K+1)/2 entries in the scaled convention documented above.
#'    diagnostics reports convergence; debiasing$status is "success",
#'    "failed", "skipped" (requested but unconverged), or "not_requested".
#'    A failed correction retains the raw fit; an invalid input or numerical
#'    fitting failure raises an error. Check status before using the results.
#'
#' The correction uses U_anchor = U_hat in the fitted coordinate frame.
#' U_hat is not corrected, and no truth-aligned confidence intervals are
#' supplied. Fitting uses no RNG, data generation, or Monte Carlo repetitions.
#' Example (A already observed): fit <- gfm_fit(A, K = 2, model = "logit")
#' Example (known sigma): gfm_fit(A, K = 2, model = "gaussian", sigma = 0.5)
gfm_fit <- function(A, K, model, sigma = NULL, control = list(), debias = TRUE) {
  if (!is.logical(debias) || length(debias) != 1L || is.na(debias)) {
    stop("debias must be TRUE or FALSE.", call. = FALSE)
  }
  A <- .gfm_validate_matrix(A)
  config <- .gfm_config(nrow(A), K, model, sigma, control)
  A <- .gfm_validate_matrix(A, config$model)
  raw <- tryCatch(.gfm_run_pgd_network(A, config), error = function(e) {
    stop("GFM fitting failed: ", conditionMessage(e), call. = FALSE)
  })
  if (raw$convergence_code == 2L) {
    stop("GFM fitting failed: non-finite objective or projected-gradient residual.",
         call. = FALSE)
  }
  labels <- rownames(A)
  coordinates <- paste0("component_", seq_len(config$k))
  dimnames(raw$U_hat) <- list(labels, coordinates)
  dimnames(raw$Lambda_hat) <- list(coordinates, coordinates)
  vnames <- .gfm_lambda_names(config$k)
  vhat <- stats::setNames(.gfm_get_lambda_v(raw$Lambda_hat), vnames)
  diagnostics <- raw[setdiff(names(raw), c("U_hat", "Lambda_hat"))]
  diagnostics$converged <- raw$convergence_code == 0L
  diagnostics$convergence_message <- switch(as.character(raw$convergence_code),
    "0" = "Projected-gradient stopping criterion satisfied.",
    "1" = "Maximum iterations reached before convergence.",
    "2" = "Non-finite objective or projected-gradient residual."
  )
  # A planning estimate, not an exact peak-RSS guarantee. Additional n^2
  # matrices, spectral initialization, R copies and LAPACK workspace cost more.
  memory_mb <- 16 * 8 * (as.double(config$n) * config$k)^2 / 1024^2
  diagnostics$estimated_debias_workspace_mb <- memory_mb
  result <- structure(list(
    U_hat = raw$U_hat, Lambda_hat = raw$Lambda_hat, lambda_v_hat = vhat,
    Lambda_debiased = NULL, lambda_v_debiased = NULL, correction = NULL,
    efficient_information = NULL, config = config, node_labels = labels,
    diagnostics = diagnostics,
    debiasing = list(status = "not_requested", message = "Debiasing was not requested.",
                     coordinate_frame = "fitted U_hat; plug-in U_anchor = U_hat")
  ), class = "gfm_fit")
  if (!diagnostics$converged) {
    if (debias) {
      result$debiasing$status <- "skipped"
      result$debiasing$message <- "Debiasing requires a converged fit."
    }
    warning(diagnostics$convergence_message, " Raw fit returned; inspect diagnostics.", call. = FALSE)
    return(result)
  }
  if (!debias) return(result)
  # Defer these checks until after estimation so failure still returns the fit.
  correction <- tryCatch({
    edges <- as.double(config$n) * (config$n - 1) / 2
    parameters <- config$n * config$k - config$k * (config$k - 1) / 2
    if (edges < parameters) stop("Too few observed edges for the model's free parameters.")
    if (memory_mb > config$max_debias_memory_mb) {
      stop(sprintf(
        "Estimated dense debiasing workspace %.1f MiB exceeds max_debias_memory_mb = %.1f.",
        memory_mb, config$max_debias_memory_mb
      ))
    }
    # U_star is unavailable: compute the same correction formulas with the
    # estimated frame as anchor. There is no oracle rotation in this API.
    ans <- .gfm_calculate_debiasing_components(
      raw$U_hat, raw$Lambda_hat, U_anchor = raw$U_hat, config = config, A = A
    )
    if (any(!is.finite(c(ans$correction_total, ans$Seff)))) {
      stop("Non-finite bias correction or efficient information.")
    }
    ans
  }, error = function(e) e)
  if (inherits(correction, "error")) {
    result$debiasing$status <- "failed"
    result$debiasing$message <- conditionMessage(correction)
    warning("Raw fit retained; debiasing failed: ", conditionMessage(correction), call. = FALSE)
    return(result)
  }
  result$lambda_v_debiased <- vhat + correction$correction_total
  result$Lambda_debiased <- .gfm_lambda_matrix(result$lambda_v_debiased, config$k)
  dimnames(result$Lambda_debiased) <- list(coordinates, coordinates)
  result$correction <- lapply(correction[c(
    "correction_likelihood", "correction_penalty", "correction_total"
  )], stats::setNames, nm = vnames)
  result$efficient_information <- correction$Seff
  dimnames(result$efficient_information) <- list(vnames, vnames)
  extra <- correction[setdiff(names(correction), c(
    "Seff", "lambda_v_hat", names(result$correction)
  ))]
  result$diagnostics <- c(result$diagnostics, extra)
  result$debiasing$status <- "success"
  result$debiasing$message <- "Additive Lambda correction computed in the fitted coordinate frame."
  result
}

#' Compute fitted edge values for the observed nodes.
#'
#' fit: converged result from gfm_fit().
#' type: "response" for conditional edge means or "link" for U Lambda U'.
#' debiased: TRUE uses Lambda_debiased and requires successful correction;
#'    FALSE uses Lambda_hat. Both keep U_hat unchanged.
#' Returns: a symmetric n-by-n matrix in the original node order, with NA
#'    on the diagonal. This supplies neither new-node predictions nor
#'    uncertainty intervals; corrected Lambda does not imply unbiased means.
#' Example: fitted_means <- gfm_predict(fit)
#' Example: corrected_means <- gfm_predict(fit, debiased = TRUE)
gfm_predict <- function(fit, type = c("response", "link"), debiased = FALSE) {
  if (!inherits(fit, "gfm_fit")) stop("fit must be returned by gfm_fit().", call. = FALSE)
  type <- match.arg(type)
  if (!is.logical(debiased) || length(debiased) != 1L || is.na(debiased)) {
    stop("debiased must be TRUE or FALSE.", call. = FALSE)
  }
  if (!isTRUE(fit$diagnostics$converged)) {
    stop("Predictions require a converged fit; inspect fit$diagnostics.", call. = FALSE)
  }
  Lambda <- if (debiased) fit$Lambda_debiased else fit$Lambda_hat
  if (is.null(Lambda)) stop("A debiased estimate is unavailable; inspect fit$debiasing.", call. = FALSE)
  theta <- fit$U_hat %*% Lambda %*% t(fit$U_hat)
  diag(theta) <- 0
  out <- if (type == "link") theta else .gfm_link_probability(theta, fit$config)
  if (any(!is.finite(out))) stop("Non-finite fitted predictions; check the model and estimates.", call. = FALSE)
  dimnames(out) <- list(fit$node_labels, fit$node_labels)
  diag(out) <- NA_real_
  out
}

#' Print a compact fit summary.
#' x: result from gfm_fit(); ... is accepted for the print generic.
#' Displays the model, dimensions, convergence and correction status, then
#' returns x invisibly. Example: print(fit).
print.gfm_fit <- function(x, ...) {
  cat(sprintf("GFM: %s, n = %d, K = %d\n", x$config$model, x$config$n, x$config$k))
  cat("Converged:", x$diagnostics$converged, "| iterations:", x$diagnostics$n_iter, "\n")
  cat("Debiasing:", x$debiasing$status, "\n")
  cat(x$debiasing$message, "\n")
  invisible(x)
}
