# Baseline inference in node-major coordinates y_i = (u_i, alpha_i).
# Loaded by gfm_baseline.R.

# Align U and Lambda to U_star; retain fitted coordinates if U_star is absent.
gfb_align <- function(fit, U_star = NULL) {
  U <- as.matrix(fit$U_hat)
  k <- ncol(U)
  if (is.null(U_star)) {
    return(list(U_q = U, alpha_q = fit$alpha_hat,
                Lambda_q = fit$Lambda_hat, Q_transpose = diag(k),
                diagnostics = list(oracle_alignment = FALSE,
                                   min_singular_value = NA_real_,
                                   condition_number = NA_real_)))
  }
  if (!identical(dim(U), dim(U_star)) || any(!is.finite(U_star))) {
    stop("U_star must be finite and have the same dimensions as U_hat.")
  }
  if (max(abs(colMeans(U_star))) > 1e-8 ||
      max(abs(crossprod(U_star) / nrow(U_star) - diag(k))) > 1e-8) {
    stop("U_star must satisfy centering and U_star' U_star = n I.")
  }
  sv <- svd(crossprod(U, U_star) / nrow(U))
  rotation <- sv$u %*% t(sv$v)
  min_sv <- min(sv$d)
  max_sv <- max(sv$d)
  if (min_sv <= 1e-10 * max(1, max_sv)) {
    stop("Oracle alignment is rank deficient or ill conditioned.")
  }
  list(
    U_q = U %*% rotation, alpha_q = fit$alpha_hat,
    Lambda_q = .gfb_symmetrize(t(rotation) %*% fit$Lambda_hat %*% rotation),
    Q_transpose = rotation,
    diagnostics = list(
      oracle_alignment = TRUE, min_singular_value = min_sv,
      condition_number = if (min_sv > .Machine$double.eps * max_sv)
        max_sv / min_sv else Inf
    )
  )
}

# G = [1'U, diag(U'U-nI)/2, offdiag(U'U)/sqrt(2),
#      offdiag(U'U_anchor-U_anchor'U)/sqrt(2)].
# Hold the anchor fixed; all alpha derivatives are zero.
.gfb_build_D_YG <- function(U, U_anchor = U) {
  if (!identical(dim(U), dim(U_anchor))) {
    stop("U and U_anchor must have identical dimensions.")
  }
  n <- nrow(U)
  k <- ncol(U)
  m <- k + 1L
  npair <- k * (k - 1L) / 2L
  out <- matrix(0, n * m, k * (k + 1L))
  for (a in seq_len(k)) {
    rows <- seq.int(a, n * m, by = m)
    out[rows, a] <- 1
    out[rows, k + a] <- U[, a]
  }
  pair <- 0L
  if (k > 1L) for (a in seq_len(k - 1L)) {
    for (b in seq.int(a + 1L, k)) {
      pair <- pair + 1L
      ra <- seq.int(a, n * m, by = m)
      rb <- seq.int(b, n * m, by = m)
      out[ra, 2L * k + pair] <- U[, b] / sqrt(2)
      out[rb, 2L * k + pair] <- U[, a] / sqrt(2)
      out[ra, 2L * k + npair + pair] <- U_anchor[, b] / sqrt(2)
      out[rb, 2L * k + npair + pair] <- -U_anchor[, a] / sqrt(2)
    }
  }
  if (any(!is.finite(out))) stop("Non-finite constraint derivatives.")
  out
}

.gfb_debias_memory_mb <- function(n, k) {
  # Estimate peak memory for dense information and eigendecompositions.
  p <- n * (k + 1L) + k * (k + 1L) / 2L
  18 * 8 * p^2 / 1024^2
}

.gfb_information <- function(U, alpha, Lambda, U_anchor = U, config,
                             A = NULL, use_expected = FALSE) {
  n <- nrow(U)
  k <- ncol(U)
  m <- k + 1L
  d <- k * (k + 1L) / 2L
  ny <- n * m
  nc <- k * (k + 1L)
  need_mb <- .gfb_debias_memory_mb(n, k)
  limit_mb <- config$max_debias_memory_mb
  if (is.null(limit_mb)) limit_mb <- 2048
  if (!is.finite(limit_mb) || limit_mb <= 0) {
    stop("max_debias_memory_mb must be a finite positive number.")
  }
  if (need_mb > limit_mb) {
    stop(sprintf("Dense baseline inference needs approximately %.0f MiB; configured limit is %.0f MiB. Use debias=FALSE for node-only inference or deliberately raise the limit.",
                 need_mb, limit_mb))
  }
  theta <- .gfb_theta(U, alpha, Lambda)
  expected <- .gfb_expected_derivative_terms(theta, config)$information
  observed <- if (use_expected) expected else
    .gfb_observed_information_matrix(theta, config, A)
  diag(expected) <- 0
  diag(observed) <- 0
  W <- cbind(U %*% Lambda, 1)
  SYY <- matrix(0, ny, ny)
  SYL <- matrix(0, ny, d)
  SLL <- matrix(0, d, d)
  expected_blocks <- lapply(seq_len(n), function(i) matrix(0, m, m))
  for (i in seq_len(n - 1L)) {
    ii <- (i - 1L) * m + seq_len(m)
    wi <- W[i, ]
    for (j in seq.int(i + 1L, n)) {
      jj <- (j - 1L) * m + seq_len(m)
      wj <- W[j, ]
      ld <- .gfb_lambda_direction_derivative(U[i, ], U[j, ])
      v <- observed[i, j]
      ev <- expected[i, j]
      SYY[ii, ii] <- SYY[ii, ii] + v * tcrossprod(wj)
      SYY[jj, jj] <- SYY[jj, jj] + v * tcrossprod(wi)
      off <- v * tcrossprod(wj, wi)
      SYY[ii, jj] <- SYY[ii, jj] + off
      SYY[jj, ii] <- SYY[jj, ii] + t(off)
      SYL[ii, ] <- SYL[ii, ] + v * tcrossprod(wj, ld)
      SYL[jj, ] <- SYL[jj, ] + v * tcrossprod(wi, ld)
      SLL <- SLL + v * tcrossprod(ld)
      expected_blocks[[i]] <- expected_blocks[[i]] + ev * tcrossprod(wj)
      expected_blocks[[j]] <- expected_blocks[[j]] + ev * tcrossprod(wi)
    }
  }
  SYY <- .gfb_symmetrize(SYY)
  SLL <- .gfb_symmetrize(SLL)
  DYG <- .gfb_build_D_YG(U, U_anchor)
  gram <- .gfb_symmetrize(crossprod(DYG))
  gram_inv <- .gfb_ordinary_inverse_pd(gram, "Baseline constraint Gram matrix")
  Pi_nui <- .gfb_symmetrize(DYG %*% gram_inv$inverse %*% t(DYG))
  PYY <- .gfb_symmetrize(diag(ny) - Pi_nui)
  idempotence <- max(abs(PYY %*% PYY - PYY))
  annihilation <- max(abs(PYY %*% DYG))
  if (idempotence > 1e-8 ||
      annihilation > 1e-8 * max(1, max(abs(DYG)))) {
    stop("Failed to construct the finite-sample baseline constraint projector.")
  }
  projected <- .gfb_symmetrize(PYY %*% SYY %*% PYY)
  projected_inv <- .gfb_safe_pinv_psd(projected, n_drop = nc)
  # Remove roundoff in dropped directions.
  projected_inv <- .gfb_symmetrize(PYY %*% projected_inv %*% PYY)
  Seff <- .gfb_symmetrize(SLL - t(SYL) %*% projected_inv %*% SYL)
  blocks <- lapply(seq_len(n), function(i) {
    ii <- (i - 1L) * m + seq_len(m)
    SYY[ii, ii, drop = FALSE]
  })
  list(
    SYY = SYY, SYL = SYL, SLY = t(SYL), SLL = SLL, Seff = Seff,
    DYG = DYG, DGY = t(DYG), DGY_DYG = gram,
    DGY_DYG_inverse = gram_inv$inverse, PYY = PYY, Pi_nui = Pi_nui,
    projected_inv = projected_inv, blocks = blocks,
    expected_blocks = expected_blocks, use_expected = use_expected,
    DYG_gram_error_from_nI = max(abs(gram / n - diag(nc))),
    DYG_gram_condition_number = gram_inv$condition_number,
    PYY_idempotence_error = idempotence,
    PYY_DYG_annihilation_error = annihilation,
    approximate_peak_memory_mb = need_mb
  )
}

.gfb_tangent_covariance <- function(info) {
  eff <- .gfb_ordinary_inverse_pd(info$Seff, "Baseline efficient information")
  R <- info$projected_inv
  coupling <- R %*% info$SYL
  OYL <- -coupling %*% eff$inverse
  OYY <- .gfb_symmetrize(R + coupling %*% eff$inverse %*% t(coupling))
  omega <- rbind(cbind(OYY, OYL), cbind(t(OYL), eff$inverse))
  list(omega = .gfb_symmetrize(omega), efficient_inverse = eff)
}

# Derivative of Lambda(v) u with respect to lambda_v. No alpha derivative.
.gfb_w_lambda_derivative <- function(u) {
  k <- length(u)
  d <- k * (k + 1L) / 2L
  out <- matrix(0, k, d)
  for (a in seq_len(k)) out[a, a] <- sqrt(2) * u[a]
  col <- k
  if (k > 1L) for (a in seq_len(k - 1L)) {
    for (b in seq.int(a + 1L, k)) {
      col <- col + 1L
      out[a, col] <- u[b]
      out[b, col] <- u[a]
    }
  }
  out
}

# Return additive constrained second-order bias using expected information.
# Probit standardization can use observed curvature; see docs/implementation_notes.md.
.gfb_bias_terms <- function(U, alpha, Lambda, info, config,
                            retain_audit_matrices = FALSE) {
  n <- nrow(U)
  k <- ncol(U)
  m <- k + 1L
  ny <- n * m
  d <- k * (k + 1L) / 2L
  p <- ny + d
  ll <- ny + seq_len(d)
  if (!isTRUE(info$use_expected) && identical(config$model, "bernoulli_probit")) {
    stop("Second-order bias requires expected information for Probit.")
  }
  theta <- .gfb_theta(U, alpha, Lambda)
  terms <- .gfb_expected_derivative_terms(theta, config)
  W <- cbind(U %*% Lambda, 1)
  covariance <- .gfb_tangent_covariance(info)
  omega <- covariance$omega
  h_ell3 <- h_score_ell2 <- h_geometry <- numeric(p)
  max_abs_em <- 0
  for (i in seq_len(n - 1L)) {
    ii <- (i - 1L) * m + seq_len(m)
    iu <- ii[seq_len(k)]
    for (j in seq.int(i + 1L, n)) {
      jj <- (j - 1L) * m + seq_len(m)
      ju <- jj[seq_len(k)]
      indices <- c(ii, jj, ll)
      g <- c(W[j, ], W[i, ],
             .gfb_lambda_direction_derivative(U[i, ], U[j, ]))
      q <- as.numeric(crossprod(g, omega[indices, indices, drop = FALSE] %*% g))
      # Hessian_theta has U_i/U_j and U_i/Lambda, U_j/Lambda blocks only.
      r <- 2 * sum(Lambda * omega[iu, ju, drop = FALSE]) +
        2 * sum(.gfb_w_lambda_derivative(U[j, ]) * omega[iu, ll, drop = FALSE]) +
        2 * sum(.gfb_w_lambda_derivative(U[i, ]) * omega[ju, ll, drop = FALSE])
      e3 <- terms$expected_ell3[i, j]
      em <- terms$expected_score_ell2[i, j]
      fisher <- terms$information[i, j]
      h_ell3[indices] <- h_ell3[indices] + 0.5 * e3 * q * g
      h_score_ell2[indices] <- h_score_ell2[indices] + em * q * g
      h_geometry[indices] <- h_geometry[indices] - 0.5 * fisher * r * g
      max_abs_em <- max(max_abs_em, abs(em))
    }
  }
  h <- h_ell3 + h_score_ell2 + h_geometry
  # Only Gram constraints have nonzero Hessians; alpha rows/columns are zero.
  constraint_trace <- numeric(k * (k + 1L))
  cov_u_sum <- matrix(0, k, k)
  for (i in seq_len(n)) {
    iu <- (i - 1L) * m + seq_len(k)
    cov_u_sum <- cov_u_sum + omega[iu, iu, drop = FALSE]
  }
  constraint_trace[k + seq_len(k)] <- diag(cov_u_sum)
  pair <- 0L
  if (k > 1L) for (a in seq_len(k - 1L)) {
    for (b in seq.int(a + 1L, k)) {
      pair <- pair + 1L
      constraint_trace[2L * k + pair] <- sqrt(2) * cov_u_sum[a, b]
    }
  }
  c_y <- as.vector(-0.5 * info$DYG %*% info$DGY_DYG_inverse %*% constraint_trace)
  c_full <- c(c_y, rep(0, d))
  I_c <- c(as.vector(info$SYY %*% c_y), as.vector(info$SLY %*% c_y))
  likelihood_bias <- as.vector(omega %*% h)
  constraint_bias <- c_full - as.vector(omega %*% I_c)
  full_bias <- likelihood_bias + constraint_bias
  correction_likelihood <- -likelihood_bias[ll]
  correction_constraint <- -constraint_bias[ll]
  full_I <- rbind(cbind(info$SYY, info$SYL), cbind(info$SLY, info$SLL))
  # Check algebraic identities, not normality.
  omega_I <- omega %*% full_I
  covariance_identity <- max(abs(omega_I %*% omega - omega)) /
    max(1, max(abs(omega)))
  constraint_identity <- max(abs(crossprod(info$DYG, full_bias[seq_len(ny)]) +
                                  0.5 * constraint_trace))
  if (covariance_identity > 1e-6 ||
      constraint_identity > 1e-7 * max(1, abs(constraint_trace))) {
    stop("Second-order covariance or constraint bias identity failed numerical validation.")
  }
  out <- list(
    method = "constrained_second_order", status = "ok",
    correction_likelihood = correction_likelihood,
    correction_constraint = correction_constraint,
    correction_total = correction_likelihood + correction_constraint,
    correction_ell3 = -as.vector(omega[ll, , drop = FALSE] %*% h_ell3),
    correction_score_ell2 = -as.vector(omega[ll, , drop = FALSE] %*% h_score_ell2),
    correction_predictor_curvature = -as.vector(omega[ll, , drop = FALSE] %*% h_geometry),
    constraint_trace = constraint_trace,
    max_abs_expected_score_ell2 = max_abs_em,
    covariance_identity_error = covariance_identity,
    constraint_bias_identity_error = constraint_identity,
    full_parameter_bias = full_bias,
    caveat = "Second-order correction derived from constrained likelihood equations; differs from the internally inconsistent literal PDF bias signs. Small-sample calibration is not guaranteed."
  )
  if (retain_audit_matrices) {
    out$omega <- omega
    out$information <- full_I
    out$h <- h
    out$c <- c_full
  }
  out
}

# Node-block bias with ell = log p and the exact finite-sample projector.
# B1 is minus the PDF's +0.5 ell''' contraction; BM is the score-Hessian term.
# Add T(B1)-T(BM), counting BM once. BP keeps its negative sign;
# its mapped additive correction enters with a plus. See docs/implementation_notes.md.
.gfb_bias_terms_pdf_bp_bm <- function(U, alpha, Lambda, info, config) {
  n <- nrow(U)
  k <- ncol(U)
  m <- k + 1L
  if (!isTRUE(info$use_expected) && identical(config$model, "bernoulli_probit")) {
    stop("BP/BM bias terms require expected information for Probit.")
  }
  terms <- .gfb_expected_derivative_terms(.gfb_theta(U, alpha, Lambda), config)
  W <- cbind(U %*% Lambda, 1)
  block_inv <- lapply(seq_len(n), function(i) {
    .gfb_ordinary_inverse_pd(info$expected_blocks[[i]],
      paste("Expected node", i, "information"))$inverse
  })
  # s[i,j] = w_i' V_j w_i in node-major y=(u,alpha) coordinates.
  s <- vapply(block_inv, function(V) rowSums((W %*% V) * W), numeric(n))
  contract <- function(coef) {
    diag(coef) <- 0
    by <- (coef * (s + t(s))) %*% W
    MU <- (coef * s) %*% U
    bl <- sqrt(2) * colSums(U * MU)
    if (k > 1L) for (a in seq_len(k - 1L)) {
      for (b in seq.int(a + 1L, k)) {
        bl <- c(bl, sum(U[, a] * MU[, b]) + sum(U[, b] * MU[, a]))
      }
    }
    list(Y = as.vector(t(by)), Lambda = bl)
  }
  B1 <- contract(-0.5 * terms$expected_ell3)
  BM <- contract(terms$expected_score_ell2)
  inverse_sum <- Reduce(`+`, block_inv)
  BP_Y <- as.vector(vapply(seq_len(n), function(i) {
    as.vector(-0.5 * (inverse_sum - block_inv[[i]]) %*% c(U[i, ], 0))
  }, numeric(m)))
  eff <- .gfb_ordinary_inverse_pd(info$Seff, "Baseline efficient information")
  lambda_block <- .gfb_ordinary_inverse_pd(info$SLL, "Lambda information block")
  T <- function(B) as.vector(eff$inverse %*%
    (B$Lambda - info$SLY %*% (info$projected_inv %*% B$Y)))
  correction_B1 <- T(B1)
  correction_BM <- -T(BM)
  Delta_Y <- info$DYG %*% (info$DGY %*% BP_Y) / n^2
  correction_BP <- as.vector(lambda_block$inverse %*% (info$SLY %*% Delta_Y))
  list(
    method = "pdf_bp_bm", status = "ok",
    correction_likelihood = correction_B1 + correction_BM,
    correction_ell3 = correction_B1, correction_B1 = correction_B1,
    correction_score_ell2 = correction_BM, correction_BM = correction_BM,
    correction_constraint = correction_BP, correction_BP = correction_BP,
    correction_total = correction_B1 + correction_BM + correction_BP,
    B1_Y = B1$Y, B1_Lambda = B1$Lambda,
    BM_Y = BM$Y, BM_Lambda = BM$Lambda, BP_Y = BP_Y,
    BP_Lambda = numeric(length(B1$Lambda)),
    max_abs_expected_score_ell2 = max(abs(terms$expected_score_ell2)),
    projector = "I-D(D'D)^-1D'",
    caveat = paste("Retains the PDF BP vector and DD' BP/n^2 structure, adds BM once,",
      "and normalizes additive signs for ell=log p. Uses the reviewed exact",
      "finite-sample projector. Increasing-dimension validity is not established by implementation checks.")
  )
}

# Compute node and global inference for a converged interior fit.
# Global debiasing requires the true U_star anchor.
gfb_infer <- function(fit, U_star = NULL, alpha_star = NULL,
                      Lambda_star = NULL, nodes = seq_len(nrow(fit$U_hat)),
                      debias = TRUE,
                      debias_method = c("pdf_bp_bm", "constrained_second_order")) {
  debias_method <- match.arg(debias_method)
  if (!inherits(fit, "gfb_fit")) stop("fit must be a gfb_fit object.")
  if (!identical(fit$diagnostics$convergence_code, 0L) &&
      !identical(fit$diagnostics$convergence_code, 0)) {
    stop("Inference requires a converged fit; inspect fit$diagnostics.")
  }
  if (length(fit$diagnostics$one_sided_nodes)) {
    stop("No interior likelihood inference is available: one-sided node observations imply a boundary alpha optimum; inspect fit$diagnostics$one_sided_nodes.")
  }
  if (!is.logical(debias) || length(debias) != 1L || is.na(debias)) {
    stop("debias must be TRUE or FALSE.")
  }
  if (debias && is.null(U_star)) {
    stop("Global debiasing requires the true U_star oracle anchor under this simulation convention.")
  }
  if (!is.null(Lambda_star) && is.null(U_star)) {
    stop("Comparing Lambda_star requires U_star to align the fitted coordinate frame.")
  }
  config <- fit$config
  u_margin <- config$M_bound_U - max(sqrt(rowSums(fit$U_hat^2)))
  a_margin <- config$M_bound_alpha - max(abs(fit$alpha_hat))
  l_margin <- config$M_bound_Lambda - max(abs(eigen(
    .gfb_symmetrize(fit$Lambda_hat), symmetric = TRUE, only.values = TRUE)$values))
  margins <- c(u_margin, a_margin, l_margin)
  bounds <- c(config$M_bound_U, config$M_bound_alpha, config$M_bound_Lambda)
  if (any(!is.finite(margins)) || any(margins <= 1e-7 * pmax(1, bounds))) {
    stop("Asymptotic inference requires interior inequality constraints; a U, alpha, or Lambda bound is active.")
  }
  aligned <- gfb_align(fit, U_star)
  U <- aligned$U_q
  alpha <- aligned$alpha_q
  Lambda <- aligned$Lambda_q
  n <- nrow(U)
  k <- ncol(U)
  m <- k + 1L
  d <- k * (k + 1L) / 2L
  if (!is.numeric(nodes) || !length(nodes) || any(!is.finite(nodes)) ||
      any(nodes != floor(nodes)) || any(nodes < 1L | nodes > n) || anyDuplicated(nodes)) {
    stop("nodes must be distinct node indices between 1 and n.")
  }
  nodes <- as.integer(nodes)
  if (!is.null(alpha_star) && (length(alpha_star) != n || any(!is.finite(alpha_star)))) {
    stop("alpha_star must be finite with one value per node.")
  }
  if (!is.null(Lambda_star) && (!identical(dim(Lambda_star), c(k, k)) ||
      any(!is.finite(Lambda_star)) || max(abs(Lambda_star - t(Lambda_star))) > 1e-10)) {
    stop("Lambda_star must be a finite symmetric K x K matrix.")
  }
  theta <- .gfb_theta(U, alpha, Lambda)
  V <- .gfb_observed_information_matrix(theta, config, fit$A)
  diag(V) <- 0
  W <- cbind(U %*% Lambda, 1)
  z_names <- c(paste0("z_u", seq_len(k)), "z_alpha_joint")
  rows <- vector("list", length(nodes))
  node_information <- vector("list", length(nodes))
  for (pos in seq_along(nodes)) {
    i <- nodes[pos]
    z <- rep(NA_real_, m)
    marginal <- chisq <- NA_real_
    info_i <- .gfb_symmetrize(crossprod(W, V[i, ] * W))
    node_information[[pos]] <- info_i
    attempt <- tryCatch({
      inv <- .gfb_ordinary_inverse_pd(info_i, paste("Node", i, "information"))
      if (!is.null(U_star) && !is.null(alpha_star)) {
        err <- c(U[i, ] - U_star[i, ], alpha[i] - alpha_star[i])
        z <- as.vector(.gfb_safe_sqrt_psd(info_i) %*% err)
        chisq <- sum(z^2)
      }
      if (!is.null(alpha_star)) {
        marginal <- (alpha[i] - alpha_star[i]) / sqrt(inv$inverse[m, m])
      }
      list(status = "ok", message = "")
    }, error = function(e) list(status = "failed", message = conditionMessage(e)))
    row <- data.frame(node = i, status = attempt$status, message = attempt$message,
                      stringsAsFactors = FALSE)
    for (a in seq_along(z_names)) row[[z_names[a]]] <- z[a]
    row$alpha_marginal_z <- marginal
    row$alpha_hat <- alpha[i]
    row$alpha_error <- if (is.null(alpha_star)) NA_real_ else alpha[i] - alpha_star[i]
    for (a in seq_len(k)) {
      row[[paste0("u_hat", a)]] <- U[i, a]
      row[[paste0("u_error", a)]] <- if (is.null(U_star)) NA_real_ else U[i, a] - U_star[i, a]
    }
    row$chi_sq_Y <- chisq
    row$df <- m
    rows[[pos]] <- row
  }
  node_table <- do.call(rbind, rows)
  names(node_information) <- as.character(nodes)
  raw_v <- .gfb_get_lambda_v(Lambda)
  true_v <- if (is.null(Lambda_star)) rep(NA_real_, d) else .gfb_get_lambda_v(Lambda_star)
  lambda <- list(status = "skipped", message = "debias=FALSE; global inference skipped.",
                 method = debias_method, lambda_v_hat = raw_v,
                 lambda_v_debiased = rep(NA_real_, d), raw_error = raw_v - true_v,
                 debiased_error = rep(NA_real_, d), z_raw = rep(NA_real_, d),
                 z_debiased = rep(NA_real_, d), chi_sq_raw = NA_real_,
                 chi_sq_debiased = NA_real_, df = d,
                 correction_likelihood = rep(NA_real_, d),
                 correction_ell3 = rep(NA_real_, d),
                 correction_score_ell2 = rep(NA_real_, d),
                 correction_BM = rep(NA_real_, d),
                 correction_constraint = rep(NA_real_, d),
                 correction_total = rep(NA_real_, d))
  diagnostic <- list(alignment = aligned$diagnostics, fit = fit$diagnostics,
                     node_standardization = "Observed block information; joint alpha coordinate mixes all y coordinates.",
                     alpha_standardization = "Marginal alpha error divided by sqrt((I_ii^-1)[K+1,K+1]).")
  if (isTRUE(debias)) {
    global <- tryCatch({
      anchor <- U_star
      observed_info <- .gfb_information(U, alpha, Lambda, anchor, config, fit$A)
      observed_inverse <- .gfb_ordinary_inverse_pd(observed_info$Seff, "Observed efficient information")
      bias_info <- if (identical(config$model, "bernoulli_probit")) {
        .gfb_information(U, alpha, Lambda, anchor, config, fit$A, use_expected = TRUE)
      } else observed_info
      correction <- if (identical(debias_method, "pdf_bp_bm")) {
        .gfb_bias_terms_pdf_bp_bm(U, alpha, Lambda, bias_info, config)
      } else .gfb_bias_terms(U, alpha, Lambda, bias_info, config)
      corrected_v <- raw_v + correction$correction_total
      half <- .gfb_safe_sqrt_psd(observed_info$Seff)
      z_raw <- as.vector(half %*% (raw_v - true_v))
      z_debiased <- as.vector(half %*% (corrected_v - true_v))
      answer <- lambda
      answer$status <- "ok"
      answer$message <- ""
      answer$lambda_v_debiased <- corrected_v
      answer$Lambda_debiased <- .gfb_lambda_matrix(corrected_v, k)
      answer$debiased_error <- corrected_v - true_v
      answer$z_raw <- z_raw
      answer$z_debiased <- z_debiased
      answer$chi_sq_raw <- if (all(is.finite(z_raw))) sum(z_raw^2) else NA_real_
      answer$chi_sq_debiased <- if (all(is.finite(z_debiased))) sum(z_debiased^2) else NA_real_
      for (nm in c("correction_likelihood", "correction_ell3", "correction_score_ell2",
                   "correction_constraint", "correction_total")) {
        answer[[nm]] <- correction[[nm]]
      }
      answer$correction_BM <- correction$correction_score_ell2
      answer$efficient_information <- observed_info$Seff
      answer$bias_diagnostics <- correction
      list(lambda = answer, diagnostics = list(
        DYG_gram_error_from_nI = observed_info$DYG_gram_error_from_nI,
        DYG_gram_condition_number = observed_info$DYG_gram_condition_number,
        PYY_idempotence_error = observed_info$PYY_idempotence_error,
        PYY_DYG_annihilation_error = observed_info$PYY_DYG_annihilation_error,
        Seff_min_eigenvalue = observed_inverse$min_eigenvalue,
        Seff_condition_number = observed_inverse$condition_number,
        approximate_peak_memory_mb = observed_info$approximate_peak_memory_mb,
        covariance_identity_error = correction$covariance_identity_error,
        constraint_bias_identity_error = correction$constraint_bias_identity_error
      ))
    }, error = function(e) {
      answer <- lambda
      answer$status <- "failed"
      answer$message <- conditionMessage(e)
      list(lambda = answer, diagnostics = list(error = conditionMessage(e)))
    })
    lambda <- global$lambda
    diagnostic$global <- global$diagnostics
  }
  failures <- sum(node_table$status != "ok") + as.integer(lambda$status == "failed")
  status <- if (!failures) "ok" else if (all(node_table$status != "ok") &&
                                          lambda$status != "ok") "failed" else "partial"
  list(status = status,
       message = if (failures) "Some inference components failed; inspect component statuses." else "",
       aligned = aligned, nodes = node_table, node_information = node_information,
       lambda = lambda, diagnostics = diagnostic)
}
