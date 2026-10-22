# Base-R likelihood, matrix and validation helpers.

.gfb_normalize_model_name <- function(model) {
    key <- tolower(gsub("[- ]", "_", model))
    aliases <- c(logit = "bernoulli_logit", bernoulli_logit = "bernoulli_logit", probit = "bernoulli_probit", 
        bernoulli_probit = "bernoulli_probit", gaussian = "gaussian", normal = "gaussian", poisson = "poisson")
    if (!key %in% names(aliases)) {
        stop("Unknown model: ", model, ". Choose bernoulli_logit, bernoulli_probit, gaussian, or poisson.")
    }
    unname(aliases[[key]])
}

.gfb_symmetrize <- function(A) (A + t(A))/2

.gfb_sigmoid <- function(x) {
    out <- numeric(length(x))
    positive <- x >= 0
    out[positive] <- 1/(1 + exp(-x[positive]))
    exp_x <- exp(x[!positive])
    out[!positive] <- exp_x/(1 + exp_x)
    dim(out) <- dim(x)
    out
}

.gfb_stable_log1pexp <- function(x) {
    pmax(x, 0) + log1p(exp(-abs(x)))
}

# Drop n_drop smallest-magnitude eigenvalues, then invert the positive remainder.
.gfb_safe_pinv_psd <- function(A, n_drop = 0L) {
    A <- .gfb_symmetrize(A)
    eg <- eigen(A, symmetric = TRUE)
    vals <- eg$values
    p <- length(vals)
    if (any(!is.finite(vals))) {
        stop("Non-finite information eigenvalue.")
    }
    if (length(n_drop) != 1L || !is.finite(n_drop) || n_drop < 0L || n_drop > p || n_drop != floor(n_drop)) {
        stop("n_drop must be an integer between 0 and the matrix dimension.")
    }
    if (n_drop == p) {
        return(matrix(0, p, p))
    }
    keep <- rep(TRUE, p)
    if (n_drop > 0L) {
        keep[order(abs(vals))[seq_len(n_drop)]] <- FALSE
    }
    if (any(vals[keep] <= 0)) {
        stop("Remaining eigenvalues must be strictly positive.")
    }
    vals_inv <- numeric(p)
    vals_inv[keep] <- 1/vals[keep]
    .gfb_symmetrize(eg$vectors %*% diag(vals_inv, p) %*% t(eg$vectors))
}

.gfb_assert_positive_definite <- function(A, label = "Information matrix", rel_tol = 1e-10) {
    values <- eigen(.gfb_symmetrize(A), symmetric = TRUE, only.values = TRUE)$values
    if (!length(values) || any(!is.finite(values))) {
        stop(label, " has non-finite or missing eigenvalues.")
    }
    scale <- max(abs(values))
    if (scale == 0 || min(values) <= rel_tol * scale) {
        stop(sprintf("%s is not numerically positive definite: min %.4e, max scale %.4e.", label, min(values), 
            scale))
    }
    invisible(values)
}

.gfb_ordinary_inverse_pd <- function(A, label = "Information matrix", rel_tol = 1e-10) {
    A <- .gfb_symmetrize(A)
    values <- .gfb_assert_positive_definite(A, label, rel_tol)
    inverse <- .gfb_symmetrize(chol2inv(chol(A)))
    inverse_residual <- max(abs(A %*% inverse - diag(nrow(A))))
    list(inverse = inverse, eigenvalues = values, min_eigenvalue = min(values), max_eigenvalue = max(values), 
        condition_number = max(values)/min(values), inverse_residual = inverse_residual)
}

# Vectorize as diag(Lambda)/sqrt(2), then upper-triangle entries by row.
.gfb_get_lambda_v <- function(Lambda) {
    k <- nrow(Lambda)
    out <- diag(Lambda)/sqrt(2)
    if (k >= 2L) {
        for (a in seq_len(k - 1L)) {
            for (b in seq.int(a + 1L, k)) out <- c(out, Lambda[a, b])
        }
    }
    as.vector(out)
}

.gfb_lambda_direction_derivative <- function(ui, uj) {
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

.gfb_proximal_projection_Lambda <- function(Lambda, bound) {
    eg <- eigen(.gfb_symmetrize(Lambda), symmetric = TRUE)
    vals <- pmax(pmin(eg$values, bound), -bound)
    .gfb_symmetrize(eg$vectors %*% diag(vals, length(vals)) %*% t(eg$vectors))
}

.gfb_probit_ratios <- function(theta) {
    log_phi <- dnorm(theta, log = TRUE)
    log_p1 <- pnorm(theta, log.p = TRUE)
    log_p0 <- pnorm(theta, lower.tail = FALSE, log.p = TRUE)
    list(p1 = pnorm(theta), p0 = pnorm(theta, lower.tail = FALSE), a = exp(log_phi - log_p1), b = exp(log_phi - 
        log_p0), log_phi = log_phi, log_p1 = log_p1, log_p0 = log_p0)
}

.gfb_link_probability <- function(theta, config) {
    switch(config$model, bernoulli_logit = .gfb_sigmoid(theta), bernoulli_probit = pnorm(theta), poisson = exp(theta), 
        gaussian = theta)
}

.gfb_link_inverse <- function(value, config, eps = config$probability_eps) {
    switch(config$model, bernoulli_logit = {
        p <- pmax(pmin(value, 1 - eps), eps)
        log(p/(1 - p))
    }, bernoulli_probit = qnorm(pmax(pmin(value, 1 - eps), eps)), poisson = log(pmax(value, 0.1)), gaussian = value)
}

.gfb_negative_loglikelihood_edges <- function(theta, y, config) {
    switch(config$model, bernoulli_logit = sum(.gfb_stable_log1pexp(theta) - y * theta), bernoulli_probit = {
        log_p1 <- pnorm(theta, log.p = TRUE)
        log_p0 <- pnorm(theta, lower.tail = FALSE, log.p = TRUE)
        -sum(ifelse(y > 0.5, log_p1, log_p0))
    }, gaussian = sum((theta - y)^2)/(2 * config$sigma^2), poisson = sum(exp(theta) - y * theta))
}

.gfb_negative_score_matrix <- function(theta, A, config) {
    switch(config$model, bernoulli_logit = .gfb_sigmoid(theta) - A, bernoulli_probit = {
        r <- .gfb_probit_ratios(theta)
        out <- ifelse(A > 0.5, -r$a, r$b)
        dim(out) <- dim(theta)
        out
    }, gaussian = (theta - A)/config$sigma^2, poisson = exp(theta) - A)
}

.gfb_expected_derivative_terms <- function(theta, config) {
    zero <- 0 * theta
    if (config$model == "bernoulli_logit") {
        p <- .gfb_sigmoid(theta)
        information <- p * (1 - p)
        expected_ell3 <- -information * (1 - 2 * p)
        return(list(mean_parameter = p, information = information, expected_ell3 = expected_ell3, expected_score_ell2 = zero, 
            observed_information_y1 = information, observed_information_y0 = information))
    }
    if (config$model == "bernoulli_probit") {
        r <- .gfb_probit_ratios(theta)
        a <- r$a
        b <- r$b
        ell2_1 <- -theta * a - a^2
        ell3_1 <- a * (theta^2 - 1) + 3 * theta * a^2 + 2 * a^3
        ell2_0 <- theta * b - b^2
        ell3_0 <- b * (1 - theta^2) + 3 * theta * b^2 - 2 * b^3
        information <- exp(2 * r$log_phi - r$log_p1 - r$log_p0)
        return(list(mean_parameter = r$p1, information = information, expected_ell3 = information * (3 * 
            theta + 2 * a - 2 * b), expected_score_ell2 = information * (b - a - theta), observed_information_y1 = -ell2_1, 
            observed_information_y0 = -ell2_0, ell3_y1 = ell3_1, ell3_y0 = ell3_0))
    }
    if (config$model == "gaussian") {
        return(list(mean_parameter = theta, information = zero + 1/config$sigma^2, expected_ell3 = zero, 
            expected_score_ell2 = zero, observed_information_y1 = zero + 1/config$sigma^2, observed_information_y0 = zero + 
                1/config$sigma^2))
    }
    intensity <- exp(theta)
    list(mean_parameter = intensity, information = intensity, expected_ell3 = -intensity, expected_score_ell2 = zero, 
        observed_information_y1 = intensity, observed_information_y0 = intensity)
}

.gfb_observed_information_matrix <- function(theta, config, A = NULL) {
    if (config$model != "bernoulli_probit") {
        return(.gfb_expected_derivative_terms(theta, config)$information)
    }
    if (is.null(A)) {
        stop("A is required for the observed probit information matrix.")
    }
    if (!identical(dim(A), dim(theta))) {
        stop("A and theta must have the same dimensions.")
    }
    r <- .gfb_probit_ratios(theta)
    information_y1 <- r$a * (theta + r$a)
    information_y0 <- r$b * (r$b - theta)
    information <- ifelse(A > 0.5, information_y1, information_y0)
    dim(information) <- dim(theta)
    pmax(information, 0)
}

.gfb_spectral_keep_indices <- function(singular_values, k, offset) {
    p <- length(singular_values)
    threshold_index <- k + offset + 1L
    threshold <- if (threshold_index <= p) {
        singular_values[threshold_index] + 1e-06
    }
    else {
        0
    }
    keep <- which(singular_values > threshold)
    if (length(keep) < k) 
        keep <- seq_len(min(k, p))
    keep
}

.gfb_validate_matrix <- function(A, model = NULL) {
    if (!is.matrix(A) || !is.numeric(A) || nrow(A) != ncol(A) || nrow(A) < 2L) {
        stop("A must be a numeric square matrix with at least two nodes.", call. = FALSE)
    }
    rn <- rownames(A)
    cn <- colnames(A)
    if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
        stop("Row and column node names must agree in the same order.", call. = FALSE)
    }
    labels <- if (!is.null(rn)) 
        rn
    else cn
    if (!is.null(labels) && (anyNA(labels) || any(!nzchar(labels)) || anyDuplicated(labels))) {
        stop("Node names must be nonempty and unique.", call. = FALSE)
    }
    diag(A) <- 0
    if (any(!is.finite(A))) 
        stop("All off-diagonal edges must be finite and observed.", call. = FALSE)
    if (!all(A == t(A))) 
        stop("A must be symmetric; edges will not be silently averaged.", call. = FALSE)
    if (!is.null(model)) {
        values <- A[upper.tri(A)]
        if (model %in% c("bernoulli_logit", "bernoulli_probit") && any(values != 0 & values != 1)) {
            stop("Logit and probit require binary (0/1) edges.", call. = FALSE)
        }
        if (model == "poisson" && any(values < 0 | values != floor(values))) {
            stop("Poisson requires nonnegative integer edges.", call. = FALSE)
        }
    }
    if (!is.null(labels)) 
        dimnames(A) <- list(labels, labels)
    A
}

# Reconstruct symmetric Lambda from lambda_v.
.gfb_lambda_matrix <- function(v, K) {
    if (length(v) != K * (K + 1L)/2L) 
        stop("Invalid lambda_v length.")
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
