# Saved-result plotting for the GFM with additive node effects.
# Sourcing defines functions only; fitting is never run.

.gfb_plot_bind <- function(x) {
  if (!length(x)) return(data.frame())
  do.call(rbind, x)
}

# Lambda_v stores diagonal / sqrt(2), followed by upper-triangle entries by row.
.gfb_plot_coordinates <- function(K) {
  out <- data.frame(row = seq_len(K), col = seq_len(K), multiplier = sqrt(2))
  if (K > 1L) for (q in seq_len(K - 1L)) for (v in seq.int(q + 1L, K)) {
    out <- rbind(out, data.frame(row = q, col = v, multiplier = 1))
  }
  out
}

.gfb_plot_required <- function(d, fields, label) {
  missing <- setdiff(fields, names(d))
  if (length(missing)) stop(label, " is missing: ", paste(missing, collapse = ", "))
}

.gfb_plot_flag <- function(d, name, allowed) {
  if (!name %in% names(d)) return(rep(TRUE, nrow(d)))
  d[[name]] %in% allowed
}

# Read per-replication results from native RDS or historical CSV files.
gfb_read_plot_data <- function(input) {
  if (length(input) != 1L || !file.exists(input)) stop("input must name results.rds or a CSV directory.")
  input <- normalizePath(input, winslash = "/", mustWork = TRUE)
  if (!dir.exists(input)) {
    sim <- readRDS(input)
    if (!is.list(sim) || !identical(sim$schema_version, "gfm-baseline-v1")) {
      stop("Expected a native gfm-baseline-v1 results.rds.")
    }
    s <- sim$settings
    t1 <- sim$node_summary; t2 <- sim$lambda_summary
    if (!is.data.frame(t1) || !is.data.frame(t2)) stop("Missing native summary data frames.")
    for (stage in c("t1", "t2")) {
      d <- get(stage)
      if (nrow(d)) {
        .gfb_plot_required(d, c("model", "replicate", "status", "fit_converged"), stage)
        d$n <- s$n; d$K <- s$K
        d$scenario <- "saved_run"
        d$lambda_values <- paste(s$lambda_values, collapse = ",")
        d$expected <- s$n_sim
        d$valid <- d$fit_converged %in% TRUE & d$status %in% "ok"
        if (stage == "t2") for (j in seq_len(s$K * (s$K + 1L) / 2L)) {
          map <- c(paste0("z_debiased_", j), paste0("lambda_v_hat_", j),
                   paste0("lambda_v_debiased_", j), paste0("raw_error_", j), paste0("debiased_error_", j))
          names(map) <- c(paste0("z_lambda_v", j, "_debiased"),
                          paste0("lambda_v_hat_v", j), paste0("lambda_v_debiased_v", j),
                          paste0("raw_error_v", j), paste0("debiased_error_v", j))
          for (target in names(map)) if (map[[target]] %in% names(d)) d[[target]] <- d[[map[[target]]]]
        }
      }
      assign(stage, d)
    }
    if (!is.null(sim$truth$Lambda_star)) {
      expected_truth <- diag(as.numeric(s$lambda_values), nrow = s$K)
      if (!isTRUE(all.equal(sim$truth$Lambda_star, expected_truth, tolerance = 1e-10, check.attributes = FALSE))) {
        stop("Saved Lambda truth disagrees with simulation settings.")
      }
    }
    inputs <- input
  } else {
    inputs <- file.path(input, c("theorem1_records.csv", "theorem2_records.csv"))
    if (!any(file.exists(inputs))) stop("Directory must contain theorem1_records.csv and/or theorem2_records.csv.")
    t1 <- if (file.exists(inputs[1])) read.csv(inputs[1], check.names = FALSE) else data.frame()
    t2 <- if (file.exists(inputs[2])) read.csv(inputs[2], check.names = FALSE) else data.frame()
    inputs <- inputs[file.exists(inputs)]
    for (stage in c("t1", "t2")) {
      d <- get(stage)
      if (nrow(d)) {
        .gfb_plot_required(d, c("model", "scenario", "n", "K", "replicate", "lambda_values",
                                "status", "converged", "inference_status"), stage)
        d$valid <- d$status %in% "success" & d$converged %in% TRUE & d$inference_status %in% "ok"
        if (stage == "t1") d$valid <- d$valid &
          .gfb_plot_flag(d, "node_inference_status", "ok") & .gfb_plot_flag(d, "inference_valid", TRUE)
        if (stage == "t2") d$valid <- d$valid &
          .gfb_plot_flag(d, "lambda_inference_status", "ok") & .gfb_plot_flag(d, "inference_valid_debiased", TRUE)
        if ("compact_fit_md5_verified" %in% names(d)) d$valid <- d$valid & d$compact_fit_md5_verified %in% TRUE
        d$expected <- NA_integer_  # Total attempts may be absent from historical exports.
      }
      assign(stage, d)
    }
  }
  all_d <- Filter(function(d) nrow(d) > 0L, list(t1, t2))
  if (!length(all_d)) stop("No saved rows found.")
  for (d in all_d) {
    if (any(!is.finite(d$K) | d$K < 1 | d$K != floor(d$K)) ||
        any(!is.finite(d$n) | d$n < 1 | d$n != floor(d$n)) ||
        any(!is.finite(d$replicate) | d$replicate < 1 | d$replicate != floor(d$replicate))) {
      stop("Invalid K, n, or replication identifiers.")
    }
    keys <- c("model", "scenario", "n", if ("node" %in% names(d)) "node", "replicate")
    if (anyDuplicated(d[keys])) stop("Duplicate replication rows: use one experiment collection per input.")
    groups <- split(d, interaction(d$model, d$scenario, d$n, drop = TRUE))
    for (g in groups) {
      for (key in intersect(c("K", "lambda_values", "alpha_seed", "design_seed", "network_seed_base"), names(g))) {
        if (length(unique(g[[key]])) != 1L) stop("Mixed settings in one group: ", key)
      }
      lambda <- suppressWarnings(as.numeric(strsplit(as.character(g$lambda_values[1]), ",", fixed = TRUE)[[1]]))
      if (length(lambda) != g$K[1] || any(!is.finite(lambda))) stop("lambda_values must contain K finite diagonal values.")
      if ("lambda_v_order" %in% names(g) && g$K[1] == 2L &&
          any(g$lambda_v_order != "Lambda11/sqrt(2);Lambda22/sqrt(2);Lambda12")) {
        stop("Unrecognized Lambda_v coordinate ordering.")
      }
    }
  }
  # Verify the stored scale before converting Lambda_v to matrix entries.
  if (nrow(t2)) for (g in split(t2, interaction(t2$model, t2$scenario, t2$n, drop = TRUE))) {
    K <- g$K[1]; m <- K * (K + 1L) / 2L
    true_v <- c(as.numeric(strsplit(as.character(g$lambda_values[1]), ",", fixed = TRUE)[[1]]) / sqrt(2), rep(0, m - K))
    for (j in seq_len(m)) for (pair in list(c("lambda_v_hat_v", "raw_error_v"), c("lambda_v_debiased_v", "debiased_error_v"))) {
      a <- paste0(pair[1], j); b <- paste0(pair[2], j)
      if (all(c(a, b) %in% names(g))) {
        finite <- g$valid & is.finite(g[[a]]) & is.finite(g[[b]])
        if (any(abs(g[[a]][finite] - g[[b]][finite] - true_v[j]) > 1e-8 * pmax(1, abs(true_v[j])))) {
          stop("Saved Lambda estimate/error disagrees with the truth: ", a)
        }
      }
    }
  }
  # Check the saved squared norm without restandardizing.
  for (stage in c("t1", "t2")) {
    d <- get(stage)
    if (!nrow(d)) next
    for (K in unique(d$K)) {
      g <- d[d$K == K, , drop = FALSE]
      zcols <- if (stage == "t1") c(paste0("z_u", seq_len(K)), "z_alpha_joint") else paste0("z_lambda_v", seq_len(K * (K + 1L) / 2L), "_debiased")
      chi <- if (stage == "t1") "chi_sq_Y" else "chi_sq_debiased"
      if (all(c(zcols, chi) %in% names(g))) {
        z <- as.matrix(g[zcols]); norm <- rowSums(z^2)
        finite <- g$valid & is.finite(norm) & is.finite(g[[chi]])
        if (any(abs(norm[finite] - g[[chi]][finite]) > 1e-8 * pmax(1, norm[finite]))) stop("Saved joint statistic is not the squared norm: ", stage)
      }
    }
  }
  # Match T1/T2 fit identities when available.
  if (nrow(t1) && nrow(t2) && all(c("task_id", "network_md5", "compact_fit_md5") %in% intersect(names(t1), names(t2)))) {
    pos <- match(t1$task_id, t2$task_id); common <- !is.na(pos)
    for (key in c("network_md5", "compact_fit_md5")) {
      if (any(t1[[key]][common] != t2[[key]][pos[common]], na.rm = TRUE)) stop("T1/T2 fit identity mismatch: ", key)
    }
  }
  list(t1 = t1, t2 = t2,
       inputs = data.frame(file = inputs, md5 = unname(tools::md5sum(inputs)), stringsAsFactors = FALSE))
}

.gfb_plot_slug <- function(x) gsub("[^A-Za-z0-9_.-]", "_", x)
.gfb_plot_model <- function(x) {
  labels <- c(gaussian = "Gaussian", poisson = "Poisson", bernoulli_logit = "Logit", bernoulli_probit = "Probit")
  if (!x %in% names(labels)) stop("Unsupported model label: ", x)
  unname(labels[x])
}

# Export a 5.5-inch vector PDF with Latin Modern labels by default.
.gfb_plot_pdf <- function(path, build_dir, draw, heading, setting, xlab, ylab,
                          font = c("latex", "serif")) {
  font <- match.arg(font)
  labels <- list()
  add <- function(x, y, text, size = 11, align = "c", angle = 0) {
    labels[[length(labels) + 1L]] <<- list(x = x, y = y, text = text, size = size, align = align, angle = angle)
  }
  stem <- tools::file_path_sans_ext(basename(path))
  background <- if (font == "latex") file.path(build_dir, paste0(stem, "_background.pdf")) else path
  grDevices::pdf(background, width = 5.5, height = 5.5, family = "Times", pointsize = 11, useDingbats = FALSE)
  closed <- FALSE
  on.exit(if (!closed) grDevices::dev.off(), add = TRUE)
  par(family = "serif", mar = c(4.5, 4.6, 3.8, 1.2), mgp = c(2.8, 0.7, 0),
      las = 1, tck = -0.018, cex.axis = 0.98, cex.lab = 1.1,
      col.axis = "#222222", col.lab = "#222222", fg = "#222222", xaxs = "i", yaxs = "i", bty = "l")
  draw(font, add)
  p <- par("plt") * 396
  if (font == "latex") {
    xt <- pretty(par("usr")[1:2], n = 4); xt <- xt[xt >= par("usr")[1] & xt <= par("usr")[2]]
    yt <- axTicks(2); yt <- yt[yt >= par("usr")[3] & yt <= par("usr")[4]]
    axis(1, at = xt, labels = FALSE, lwd = 0, lwd.ticks = .7)
    axis(2, at = yt, labels = FALSE, lwd = 0, lwd.ticks = .7)
    for (i in seq_along(xt)) add(grconvertX(xt[i], "user", "ndc") * 396, p[3] - 15,
                                paste0("$", format(xt[i], trim = TRUE, scientific = FALSE), "$"), 10.5)
    for (i in seq_along(yt)) add(p[1] - 8, grconvertY(yt[i], "user", "ndc") * 396,
                                paste0("$", format(yt[i], trim = TRUE, scientific = FALSE), "$"), 10.5, "r")
    add(mean(p[1:2]), 377, heading, 12.1)
    add(mean(p[1:2]), 360, setting, 10.1)
    add(mean(p[1:2]), 14, xlab, 12.1)
    add(16, mean(p[3:4]), ylab, 12.1, angle = 90)
  } else {
    axis(1); axis(2)
    plain <- function(x) {
      x <- gsub("\\Lambda", "Lambda", x, fixed = TRUE)
      x <- gsub("\\alpha", "alpha", x, fixed = TRUE)
      x <- gsub("\\star", "*", x, fixed = TRUE)
      x <- gsub("\\chi", "chi", x, fixed = TRUE)
      x <- gsub("\\quad", "  ", x, fixed = TRUE)
      gsub("[{}$]", "", gsub("\\\\[A-Za-z]+", "", x))
    }
    title(main = plain(heading), cex.main = .95, font.main = 1, line = 2)
    mtext(plain(setting), side = 3, line = .8, cex = .72)
    title(xlab = plain(xlab), ylab = plain(ylab))
  }
  box(bty = "l", lwd = .7, col = "#444444")
  grDevices::dev.off(); closed <- TRUE
  if (font == "latex") {
    commands <- vapply(labels, function(a) {
      content <- sprintf("{\\fontsize{%.3fbp}{%.3fbp}\\selectfont %s}", a$size, a$size * 1.2, a$text)
      if (a$angle != 0) content <- sprintf("\\rotatebox{%d}{%s}", a$angle, content)
      sprintf("\\put(%.5f,%.5f){\\makebox(0,0)[%s]{%s}}", a$x, a$y, a$align, content)
    }, character(1))
    tex <- file.path(build_dir, paste0(stem, ".tex"))
    writeLines(c("\\documentclass{article}",
      "\\usepackage[paperwidth=5.5in,paperheight=5.5in,margin=0pt]{geometry}",
      "\\usepackage[T1]{fontenc}", "\\usepackage{lmodern,amsmath,amssymb,graphicx,xcolor}",
      "\\definecolor{ink}{HTML}{222222}", "\\pagestyle{empty}", "\\setlength{\\parindent}{0pt}",
      "\\setlength{\\topskip}{0pt}", "\\setlength{\\unitlength}{1bp}", "\\begin{document}",
      "\\color{ink}", "\\noindent\\begin{picture}(396,396)",
      paste0("\\put(0,0){\\makebox(0,0)[lb]{\\includegraphics[width=396bp,height=396bp]{", basename(background), "}}}"),
      commands, "\\end{picture}", "\\end{document}"), tex)
    old <- setwd(build_dir); on.exit(setwd(old), add = TRUE)
    log <- system2(Sys.which("pdflatex"), c("-interaction=nonstopmode", "-halt-on-error", "-no-shell-escape", shQuote(basename(tex))), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(log, "status")) && attr(log, "status") != 0) stop("LaTeX figure failed:\n", paste(tail(log, 20), collapse = "\n"))
    if (!file.copy(paste0(stem, ".pdf"), path, overwrite = FALSE)) stop("Cannot save figure: ", path)
  }
  invisible(path)
}

# Export QQ PDFs and paired histograms from saved results.
# Histograms use matrix entries; QQ plots use standardized statistics.
gfb_plot_saved <- function(input, output_dir, nodes = 1L, models = NULL, sample_sizes = NULL,
                           kinds = c("qq", "hist_estimates", "hist_errors"),
                           font = c("latex", "serif"), include_alpha_marginal = TRUE,
                           histogram_breaks = "FD") {
  font <- match.arg(font)
  if (!length(kinds) || any(!kinds %in% c("qq", "hist_estimates", "hist_errors"))) stop("Unknown plot kind.")
  kinds <- unique(kinds)
  if (font == "latex" && !nzchar(Sys.which("pdflatex"))) stop("LaTeX figures require pdflatex and lmodern; use font='serif' for base-R text.")
  if (!is.null(models)) {
    aliases <- c(logit = "bernoulli_logit", probit = "bernoulli_probit", normal = "gaussian")
    models <- tolower(models)
    hit <- models %in% names(aliases); models[hit] <- unname(aliases[models[hit]])
    if (any(!models %in% c("gaussian", "poisson", "bernoulli_logit", "bernoulli_probit"))) stop("Unknown model filter.")
  }
  saved <- gfb_read_plot_data(input)
  t1 <- saved$t1; t2 <- saved$t2
  select <- function(d) {
    if (!nrow(d)) return(d)
    if (!is.null(models)) d <- d[d$model %in% models, , drop = FALSE]
    if (!is.null(sample_sizes)) d <- d[d$n %in% sample_sizes, , drop = FALSE]
    d
  }
  t1 <- select(t1); t2 <- select(t2)
  if (nrow(t1)) {
    .gfb_plot_required(t1, "node", "T1")
    if (!identical(nodes, "all")) {
      nodes <- if (is.character(nodes)) suppressWarnings(as.numeric(strsplit(nodes, ",", fixed = TRUE)[[1]])) else nodes
      if (!length(nodes) || anyNA(nodes) || any(nodes != floor(nodes)) || any(!nodes %in% t1$node)) stop("nodes must select stored node indices, or 'all'.")
      t1 <- t1[t1$node %in% nodes, , drop = FALSE]
    }
  }
  if (!nrow(t1) && !nrow(t2)) stop("No rows match the requested filters.")
  hist_requested <- any(c("hist_estimates", "hist_errors") %in% kinds)
  if (hist_requested) {
    if (!nrow(t2)) stop("Histograms require a T2 table with saved raw and debiased Lambda estimates.")
    for (K in unique(t2$K)) {
      m <- K * (K + 1L) / 2L
      .gfb_plot_required(t2, c(paste0("lambda_v_hat_v", seq_len(m)), paste0("lambda_v_debiased_v", seq_len(m))), "Histogram input")
    }
    if (!any(t2$valid)) stop("No successful debiased Lambda records are available for paired histograms.")
  }
  if (length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) stop("Choose a new output_dir.")
  if (file.exists(output_dir)) stop("Output already exists; choose a new directory: ", output_dir)
  if (!dir.create(output_dir, recursive = TRUE)) stop("Cannot create output directory.")
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  build_dir <- file.path(output_dir, "_build"); dir.create(build_dir)
  audits <- points <- hist_data <- list()
  groups <- function(d, node = FALSE) {
    if (!nrow(d)) return(list())
    split(d, interaction(d$model, d$scenario, d$n, if (node) d$node else 1L, drop = TRUE))
  }
  setting <- function(d, node = FALSE) {
    lambda <- as.numeric(strsplit(as.character(d$lambda_values[1]), ",", fixed = TRUE)[[1]])
    paste0(.gfb_plot_model(d$model[1]), " model,\\quad ",
           if (node) paste0("node ", d$node[1], ",\\quad ") else "",
           "$\\Lambda^{\\star}=\\operatorname{diag}(", paste(lambda, collapse = ",\\,"), ")$")
  }
  audit <- function(d, kind, column, used, figure, node = NA_integer_) {
    audits[[length(audits) + 1L]] <<- data.frame(model = d$model[1], scenario = d$scenario[1],
      n = d$n[1], K = d$K[1], node = node, kind = kind, statistic = column,
      stored_rows = nrow(d), retained_rows = sum(used), excluded_rows = sum(!used),
      expected_repetitions = d$expected[1],
      missing_rows = if (is.na(d$expected[1])) NA_integer_ else max(d$expected[1] - nrow(d), 0L), figure = figure, stringsAsFactors = FALSE)
  }
  stem <- function(d, type, column, node = FALSE) {
    paste(.gfb_plot_slug(d$model[1]), .gfb_plot_slug(d$scenario[1]), paste0("n", d$n[1]),
          if (node) paste0("node", d$node[1]) else "global", type, column, sep = "_")
  }
  qq <- function(d, column, statistic, type, df = NULL, node = FALSE) {
    if (!column %in% names(d)) return(invisible(NULL))
    if (!is.numeric(d[[column]])) stop("Non-numeric statistic: ", column)
    used <- d$valid & is.finite(d[[column]])
    name <- paste0(stem(d, type, column, node), ".pdf")
    if (sum(used) < 2L) { audit(d, type, column, used, NA_character_, if(node) d$node[1] else NA_integer_); return(invisible(NULL)) }
    z <- d[used, , drop = FALSE]; ord <- order(z[[column]], z$replicate)
    probability <- (seq_len(nrow(z)) - .5) / nrow(z)
    x <- if (is.null(df)) qnorm(probability) else qchisq(probability, df)
    y <- z[[column]][ord]
    limits <- if (is.null(df)) {a <- ceiling(max(abs(c(x, y))) * 1.04 * 2) / 2; c(-a, a)} else c(0, ceiling(max(c(x, y)) * 1.04))
    draw <- function(font, add) {
      plot.new(); plot.window(limits, limits, asp = 1)
      abline(0, 1, lty = 2, lwd = .85, col = "#777777")
      points(x, y, pch = 1, cex = .48, lwd = .65, col = "#285D89")
    }
    .gfb_plot_pdf(file.path(output_dir, name), build_dir, draw,
      paste0(if (is.null(df)) "Normal Q-Q plot of $" else "Chi-square Q-Q plot of $", statistic, "$"),
      setting(d, node), if (is.null(df)) "Standard normal quantiles" else paste0("$\\chi^2_{", df, "}$ quantiles"),
      "Empirical quantiles", font)
    points[[length(points) + 1L]] <<- data.frame(figure = name, replicate = z$replicate[ord],
      rank = seq_along(x), probability = probability, theoretical = x, empirical = y)
    audit(d, type, column, used, name, if(node) d$node[1] else NA_integer_)
  }
  if ("qq" %in% kinds) {
    for (d in groups(t1, TRUE)) {
      K <- d$K[1]
      for (j in seq_len(K)) qq(d, paste0("z_u", j), paste0("Z_{Y,", j, "}"), "theorem1", node = TRUE)
      qq(d, "z_alpha_joint", paste0("Z_{Y,", K + 1L, "}"), "theorem1", node = TRUE)
      if (isTRUE(include_alpha_marginal)) qq(d, "alpha_marginal_z", "Z_{\\alpha}^{\\mathrm{marg}}", "alpha_marginal", node = TRUE)
      qq(d, "chi_sq_Y", "\\|Z_Y\\|^2", "theorem1", df = K + 1L, node = TRUE)
    }
    for (d in groups(t2)) {
      m <- d$K[1] * (d$K[1] + 1L) / 2L
      for (j in seq_len(m)) qq(d, paste0("z_lambda_v", j, "_debiased"), paste0("Z_{\\Lambda,", j, "}^{\\mathrm{db}}"), "theorem2_debiased")
      qq(d, "chi_sq_debiased", "\\|Z_{\\Lambda}^{\\mathrm{db}}\\|^2", "theorem2_debiased", df = m)
    }
  }
  if (hist_requested) for (d in groups(t2)) {
    coordinates <- .gfb_plot_coordinates(d$K[1])
    truth_diag <- as.numeric(strsplit(as.character(d$lambda_values[1]), ",", fixed = TRUE)[[1]])
    for (j in seq_len(nrow(coordinates))) {
      row_index <- coordinates$row[j]; col_index <- coordinates$col[j]; multiplier <- coordinates$multiplier[j]
      raw_col <- paste0("lambda_v_hat_v", j); db_col <- paste0("lambda_v_debiased_v", j)
      if (!is.numeric(d[[raw_col]]) || !is.numeric(d[[db_col]])) stop("Histogram estimates must be numeric.")
      used <- d$valid & is.finite(d[[raw_col]]) & is.finite(d[[db_col]])
      true <- if (row_index == col_index) truth_diag[row_index] else 0
      raw <- multiplier * d[[raw_col]][used]; debiased <- multiplier * d[[db_col]][used]
      # Share error-scale bins; shift by truth for the estimate view.
      error_breaks <- if (sum(used) >= 2L) hist(c(raw - true, debiased - true), breaks = histogram_breaks, plot = FALSE)$breaks else NULL
      ha <- if (sum(used) >= 2L) hist(raw - true, breaks = error_breaks, plot = FALSE) else NULL
      hb <- if (sum(used) >= 2L) hist(debiased - true, breaks = error_breaks, plot = FALSE) else NULL
      for (kind in intersect(c("hist_estimates", "hist_errors"), kinds)) {
        column <- paste0("Lambda", row_index, "_", col_index); name <- paste0(stem(d, kind, column), ".pdf")
        if (sum(used) < 2L) {audit(d, kind, column, used, NA_character_); next}
        shift <- if (kind == "hist_errors") true else 0
        a <- raw - shift; b <- debiased - shift; ref <- true - shift
        # Reuse both histograms after shifting the bins.
        breaks <- error_breaks + ref
        math <- paste0("\\Lambda_{", row_index, ",", col_index, "}")
        legend_tex <- c(paste0("$\\widehat{\\Lambda}_{", row_index, ",", col_index, "}", if (kind == "hist_errors") paste0("-", math, "^{\\star}") else "", "$"),
                        paste0("$\\widetilde{\\Lambda}_{", row_index, ",", col_index, "}", if (kind == "hist_errors") paste0("-", math, "^{\\star}") else "", "$"))
        draw <- function(font, add) {
          plot.new(); plot.window(range(c(breaks, ref)), c(0, max(ha$density, hb$density) * 1.22))
          rect(head(breaks, -1), 0, tail(breaks, -1), ha$density, col = adjustcolor("#E86D5A", .45), border = "#AB493B", lwd = .55)
          rect(head(breaks, -1), 0, tail(breaks, -1), hb$density, col = adjustcolor("#4C78C4", .45), border = "#285D89", lwd = .55)
          abline(v = ref, lty = 3, lwd = .8, col = "#555555")
          p <- par("plt") * 396
          for (k in 1:2) {
            x <- p[1] + 10 + (k - 1) * 135; y <- p[4] - 12
            xx <- grconvertX(c(x, x + 11) / 396, "ndc", "user")
            yy <- grconvertY(c(y - 4, y + 4) / 396, "ndc", "user")
            rect(xx[1], yy[1], xx[2], yy[2], col = adjustcolor(c("#E86D5A", "#4C78C4")[k], .45), border = c("#AB493B", "#285D89")[k], lwd = .55)
            if (font == "latex") add(x + 16, y, legend_tex[k], 10.5, "l") else text(grconvertX((x + 16) / 396, "ndc", "user"), mean(yy), c("Raw", "Debiased")[k], adj = c(0, .5), cex = .8)
          }
        }
        .gfb_plot_pdf(file.path(output_dir, name), build_dir, draw,
          paste0(if (kind == "hist_errors") "Estimation errors for $" else "Estimates of $", math, "$"),
          setting(d), if (kind == "hist_errors") "Estimation error" else "Estimate", "Density", font)
        hist_data[[length(hist_data) + 1L]] <- data.frame(figure = name, replicate = d$replicate[used],
          row = row_index, col = col_index, truth = true, raw = a, debiased = b, scale = if(kind == "hist_errors") "matrix_entry_error" else "matrix_entry")
        audit(d, kind, column, used, name)
      }
    }
  }
  result <- .gfb_plot_bind(audits)
  write.csv(result, file.path(output_dir, "figure_audit.csv"), row.names = FALSE)
  if (length(points)) write.csv(.gfb_plot_bind(points), file.path(output_dir, "qq_points.csv"), row.names = FALSE)
  if (length(hist_data)) write.csv(.gfb_plot_bind(hist_data), file.path(output_dir, "histogram_values.csv"), row.names = FALSE)
  if (!identical(saved$inputs$md5, unname(tools::md5sum(saved$inputs$file)))) stop("Input changed during plotting.")
  write.csv(saved$inputs, file.path(output_dir, "input_hashes.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
  saveRDS(list(kinds = kinds, font = font, nodes = nodes, models = models, sample_sizes = sample_sizes,
               include_alpha_marginal = include_alpha_marginal, histogram_breaks = histogram_breaks,
               qq_positions = "(rank - 0.5) / retained_rows", histogram_scale = "matrix entries",
               diagonal_multiplier = sqrt(2), reference_bands = FALSE), file.path(output_dir, "plot_settings.rds"))
  invisible(list(output_dir = output_dir, figures = result))
}
