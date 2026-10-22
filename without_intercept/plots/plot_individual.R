# Plot saved no-intercept results; sourcing only defines functions.

.gfm_plot_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < 1 || x != floor(x)) stop(name, " must be a positive integer.")
  as.integer(x)
}

# Read result/config pairs or a plotting manifest; record any excluded rows.
gfm_read_figure_cases <- function(data_dir, n = 1000,
                                  models = c("gaussian", "poisson", "logit", "probit"),
                                  lambda_settings = NULL,
                                  figure_types = c("qq", "hist_estimates", "hist_errors"),
                                  require_complete = TRUE) {
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  n <- .gfm_plot_integer(n, "n")
  figure_types <- unique(as.character(unlist(figure_types)))
  if (!length(figure_types) || !all(figure_types %in% c("qq", "hist_estimates", "hist_errors")))
    stop("Unknown figure_types.")
  if (!is.logical(require_complete) || length(require_complete) != 1L || is.na(require_complete))
    stop("require_complete must be TRUE or FALSE.")
  has_qq <- "qq" %in% figure_types
  has_hist <- any(c("hist_estimates", "hist_errors") %in% figure_types)
  models <- unique(sub("^bernoulli_", "", as.character(unlist(models))))
  supported <- c("gaussian", "poisson", "logit", "probit")
  if (!length(models) || !all(models %in% supported)) stop("Unsupported models.")
  if (!is.null(lambda_settings)) {
    if (!is.list(lambda_settings)) lambda_settings <- list(lambda_settings)
    lambda_settings <- lapply(lambda_settings, function(v) {
      v <- as.numeric(unlist(v))
      if (!length(v) || any(!is.finite(v))) stop("Invalid lambda_settings.")
      v
    })
  }
  resolve <- function(path) {
    if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) return(path)
    candidates <- c(file.path(data_dir, path), file.path(data_dir, sub("^data/", "", path)))
    present <- unique(candidates[file.exists(candidates)])
    if (!length(present)) stop("Missing manifest input: ", path)
    normalizePath(present[1])
  }
  manifest_path <- file.path(data_dir, "manifest.csv")
  if (file.exists(manifest_path)) {
    manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
    if (!all(c("csv_path", "config_path") %in% names(manifest)))
      stop("manifest.csv must identify csv_path and config_path.")
    pairs <- data.frame(csv = vapply(manifest$csv_path, resolve, ""),
                        config = vapply(manifest$config_path, resolve, ""))
  } else {
    csv <- sort(list.files(data_dir, "_results[.]csv$", recursive = TRUE, full.names = TRUE))
    pairs <- data.frame(csv = csv, config = sub("_results[.]csv$", "_config.rds", csv))
  }
  if (!nrow(pairs)) stop("No saved result/config pairs found in data_dir.")
  number_cols <- function(d, stem) {
    z <- grep(paste0("^", stem, "[0-9]+$"), names(d), value = TRUE)
    z[order(as.integer(sub(stem, "", z, fixed = TRUE)))]
  }
  cases <- list()
  for (i in seq_len(nrow(pairs))) {
    cfg <- readRDS(pairs$config[i])
    fields <- c("model", "n", "k", "n_sim", "lambda_vec", "target_node",
                "design_seed", "replication_seed_base", "truth_design")
    if (!is.list(cfg) || !all(fields %in% names(cfg))) stop("Incomplete config: ", pairs$config[i])
    cfg$model <- sub("^bernoulli_", "", cfg$model)
    if (length(cfg$model) != 1L || is.na(cfg$model) || !cfg$model %in% supported)
      stop("Invalid model in config.")
    for (field in c("n", "k", "n_sim", "target_node"))
      cfg[[field]] <- .gfm_plot_integer(cfg[[field]], paste("config", field))
    if (!is.numeric(cfg$lambda_vec) || length(cfg$lambda_vec) != cfg$k ||
        any(!is.finite(cfg$lambda_vec)) || cfg$target_node > cfg$n) stop("Invalid config dimensions.")
    if (cfg$n != n || !cfg$model %in% models) next
    if (!is.null(lambda_settings) && !any(vapply(lambda_settings,
        function(v) length(v) == cfg$k && isTRUE(all.equal(v, cfg$lambda_vec, tolerance = 1e-12)), logical(1)))) next
    K <- cfg$k; P <- as.integer(K * (K + 1L) / 2L)
    d <- read.csv(pairs$csv[i], check.names = FALSE)
    stems <- c(if (has_qq) c("stat_U_", "debiased_lambda_"),
               if (has_hist) c("lambda_raw_error_", "lambda_debiased_error_"))
    for (j in seq_along(stems)) {
      count <- if (stems[j] == "stat_U_") K else P
      if (!identical(number_cols(d, stems[j]), paste0(stems[j], seq_len(count))))
        stop("Component columns disagree with K: ", pairs$csv[i])
    }
    required <- c(unlist(lapply(stems, function(z) number_cols(d, z))),
                  if (has_qq) c("chi_sq_U", "chi_sq_debiased_lambda"))
    if (!all(c(required, "convergence_code") %in% names(d)) || !nrow(d))
      stop("Missing required columns or empty input: ", pairs$csv[i])
    for (column in required) if (is.logical(d[[column]]) && all(is.na(d[[column]])))
      d[[column]] <- as.numeric(d[[column]])
    if (!all(vapply(d[required], is.numeric, logical(1))) || !is.numeric(d$convergence_code))
      stop("Nonnumeric statistics: ", pairs$csv[i])
    converged <- !is.na(d$convergence_code) & d$convergence_code == 0L
    finite <- rowSums(!is.finite(as.matrix(d[required]))) == 0L
    keep <- converged & finite
    selection <- data.frame(row = seq_len(nrow(d)), converged = converged, finite = finite, included = keep)
    attempted <- nrow(d)
    if (require_complete && (attempted != cfg$n_sim || !all(keep)))
      stop("Require complete successful replications with finite selected statistics: ", pairs$csv[i],
           ". Set require_complete=FALSE to retain valid rows with an exclusion record.")
    if (sum(keep) < 2L) stop("Fewer than two usable replications: ", pairs$csv[i])
    if (any(!keep) || attempted != cfg$n_sim)
      warning("Selected ", sum(keep), " of ", attempted, " saved rows (planned ", cfg$n_sim,
              "). See records/row_selection.csv: ", pairs$csv[i])
    d <- d[keep, , drop = FALSE]
    if (has_qq) for (spec in list(c("chi_sq_U", "stat_U_"), c("chi_sq_debiased_lambda", "debiased_lambda_"))) {
      norm_sq <- rowSums(d[number_cols(d, spec[2])]^2)
      if (any(abs(d[[spec[1]]] - norm_sq) > 1e-8 * pmax(1, norm_sq)))
        stop("Joint statistic disagrees with components: ", pairs$csv[i])
    }
    ij <- cbind(a = seq_len(K), b = seq_len(K))
    if (K > 1L) for (a in seq_len(K - 1L)) ij <- rbind(ij, cbind(a, seq.int(a + 1L, K)))
    coordinates <- data.frame(component = seq_len(P), a = ij[, 1], b = ij[, 2],
      scale = ifelse(ij[, 1] == ij[, 2], sqrt(2), 1))
    truth <- diag(cfg$lambda_vec, nrow = K, ncol = K)
    true_entries <- truth[ij]
    # lambda_v stores diagonal entries divided by sqrt(2); restore matrix units.
    raw_error <- if (has_hist) sweep(as.matrix(d[paste0("lambda_raw_error_", seq_len(P))]), 2, coordinates$scale, "*") else NULL
    debiased_error <- if (has_hist) sweep(as.matrix(d[paste0("lambda_debiased_error_", seq_len(P))]), 2, coordinates$scale, "*") else NULL
    tag <- if (K == 2L && isTRUE(all.equal(cfg$lambda_vec, c(.5, .5)))) "equal" else
      if (K == 2L && isTRUE(all.equal(cfg$lambda_vec, c(.7, -.5)))) "mixed" else
      paste0("lambda_", paste(gsub("-", "m", gsub("[.]", "p", format(cfg$lambda_vec, trim = TRUE, scientific = FALSE, digits = 15))), collapse = "_"))
    cases[[length(cases) + 1L]] <- list(model = cfg$model, scenario = tag, cfg = cfg, d = d,
      K = K, P = P, coordinates = coordinates, true_entries = true_entries,
      raw_error = raw_error, debiased_error = debiased_error,
      raw_estimate = if (has_hist) sweep(raw_error, 2, true_entries, "+") else NULL,
      debiased_estimate = if (has_hist) sweep(debiased_error, 2, true_entries, "+") else NULL,
      attempted = attempted, selection = selection,
      csv = normalizePath(pairs$csv[i]), config_path = normalizePath(pairs$config[i]))
  }
  if (!length(cases)) stop("No saved cases match n, models, and lambda_settings.")
  keys <- vapply(cases, function(x) paste(x$model, x$cfg$n, x$K,
    paste(sprintf("%.17g", x$cfg$lambda_vec), collapse = ","), sep = "|"), "")
  if (anyDuplicated(keys)) stop("Duplicate settings: choose one coherent experiment collection.")
  design <- vapply(cases, function(x) paste(x$cfg$design_seed, x$cfg$replication_seed_base,
    x$cfg$target_node, x$cfg$truth_design, sep = "|"), "")
  if (length(unique(design)) != 1L) stop("Inputs have incompatible seeds, target nodes, or designs.")
  cases
}

# Export LaTeX QQ plots and estimate/error histograms as separate vector PDFs.
gfm_plot_individual <- function(data_dir, output_dir, n = 1000,
                                models = c("gaussian", "poisson", "logit", "probit"),
                                lambda_settings = NULL, histogram_bins = 30,
                                figure_types = c("qq", "hist_estimates", "hist_errors"),
                                overwrite = FALSE, previews = TRUE, require_complete = TRUE) {
  positive_integer <- .gfm_plot_integer
  n_bins <- positive_integer(histogram_bins, "histogram_bins")
  n <- positive_integer(n, "n")
  for (value in list(overwrite, previews)) if (!is.logical(value) || length(value) != 1L || is.na(value))
    stop("overwrite and previews must be TRUE or FALSE.")
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir))
    stop("output_dir must be a path.")
  figure_types <- unique(as.character(unlist(figure_types)))
  if (!length(figure_types) || !all(figure_types %in% c("qq", "hist_estimates", "hist_errors")))
    stop("Unknown figure_types.")
  if (!nzchar(Sys.which("pdflatex"))) stop("Install TeX Live/TinyTeX with pdflatex and lmodern.")
  if (!capabilities("cairo")) stop("R with Cairo support is required.")
  cases <- gfm_read_figure_cases(data_dir, n, models, lambda_settings, figure_types, require_complete)
  model_labels <- c(gaussian = "Gaussian", poisson = "Poisson", logit = "Logit", probit = "Probit")
  params <- list(n = n, models = models, lambda_settings = lambda_settings,
    histogram_bins = histogram_bins, figure_types = figure_types, overwrite = overwrite, previews = previews,
    require_complete = require_complete)

  # Use a fresh output folder unless overwrite is requested.
  base_dir <- output_dir
  output_dir <- base_dir
  if (!overwrite && dir.exists(base_dir) &&
      length(list.files(base_dir, all.files = TRUE, no.. = TRUE))) {
    i <- 1L
    repeat {
      output_dir <- file.path(base_dir, sprintf("run_%03d", i))
      if (!file.exists(output_dir)) break
      i <- i + 1L
    }
  }
  for (folder in c("qq", "hist_estimates", "hist_errors", "previews", "records", "latex"))
    dir.create(file.path(output_dir, folder), recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir)
  input_records <- do.call(rbind, lapply(cases, function(x) data.frame(
    model = x$model, scenario = x$scenario, n = x$cfg$n, K = x$K,
    target_node = x$cfg$target_node, planned = x$cfg$n_sim, attempted = x$attempted,
    replications = nrow(x$d), excluded = x$attempted - nrow(x$d),
    lambda = paste(x$cfg$lambda_vec, collapse = ","), csv = x$csv,
    csv_md5 = unname(tools::md5sum(x$csv)), config = x$config_path,
    config_md5 = unname(tools::md5sum(x$config_path)))))
  write.csv(input_records, file.path(output_dir, "records", "inputs.csv"), row.names = FALSE)
  selection <- do.call(rbind, lapply(cases, function(x) cbind(csv = x$csv, x$selection)))
  write.csv(selection, file.path(output_dir, "records", "row_selection.csv"), row.names = FALSE)

  build_dir <- file.path(output_dir, "latex")


  engine <- Sys.which("pdflatex")
  if (!nzchar(engine)) stop("Install TinyTeX/TeX Live with pdflatex and the lmodern package.")
  if (!capabilities("cairo")) stop("This notebook requires R with Cairo support.")
  preview_engine <- if (previews) Sys.which("pdftoppm") else ""
  if (previews && !nzchar(preview_engine)) message("pdftoppm not found; exporting PDFs without PNG previews.")
  raw_fill <- adjustcolor("#E86D5A", alpha.f=0.56)
  debiased_fill <- adjustcolor("#4C78C4", alpha.f=0.56)
  figure_records <- list()
  bin_records <- list()
  summary_records <- list()
  labels <- list()
  plot_style <- function() {
    par(family="serif", mar=c(4.5,4.6,3.8,1.2),
        mgp=c(2.8,0.7,0), las=1, tck=-0.018,
        cex.axis=0.98, cex.lab=1.10,
        col.axis="#222222", col.lab="#222222", fg="#222222",
        xaxs="i", yaxs="i", bty="l")
  }
  labels <- list()
  label_at <- function(x, y, text, size=11, align="c", angle=0) {
    labels[[length(labels)+1L]] <<- list(x=x, y=y, text=text,
      size=size, align=align, angle=angle)
  }
  xbp <- function(v) grconvertX(v, from="user", to="ndc") * 396
  ybp <- function(v) grconvertY(v, from="user", to="ndc") * 396

  # Typeset ticks in LaTeX.
  draw_axes <- function() {
    xt <- pretty(par("usr")[1:2], n=4)
    xt <- xt[xt >= par("usr")[1] & xt <= par("usr")[2]]
    yt <- axTicks(2)
    yt <- yt[yt >= par("usr")[3] & yt <= par("usr")[4]]
    axis(1, at=xt, labels=FALSE, lwd=0, lwd.ticks=0.7)
    axis(2, at=yt, labels=FALSE, lwd=0, lwd.ticks=0.7)
    box(bty="l", lwd=0.7, col="#444444")
    p <- par("plt") * 396
    xt_labels <- format(xt, trim=TRUE, scientific=FALSE)
    yt_labels <- format(yt, trim=TRUE, scientific=FALSE)
    for (i in seq_along(xt)) label_at(xbp(xt[i]), p[3]-15,
      paste0("$", xt_labels[i], "$"), size=10.78)
    for (i in seq_along(yt)) label_at(p[1]-8, ybp(yt[i]),
      paste0("$", yt_labels[i], "$"), size=10.78, align="r")
  }

  add_labels <- function(x, heading, xlab, ylab, node=FALSE) {
    p <- par("plt") * 396
    center <- mean(p[1:2])
    label_at(center, 377, heading, size=12.1)
    setting <- paste0(model_labels[[x$model]], " model,\\quad ",
      if (node) paste0("node ", x$cfg$target_node, ",\\quad ") else "",
      "$\\Lambda^{\\star}=\\operatorname{diag}(",
      paste(x$cfg$lambda_vec, collapse=",\\,"), ")$")
    label_at(center, 360, setting, size=10.1)
    label_at(center, 14, xlab, size=12.1)
    label_at(16, mean(p[3:4]), ylab, size=12.1, angle=90)
  }
  theoretical <- function(R, df = NULL) {
    p <- (seq_len(R) - 0.5) / R
    if (is.null(df)) qnorm(p) else qchisq(p, df = df)
  }
  # Match QQ axes across models/settings for the same dimension and statistic.
  qq_limits <- function(x, column, df = NULL) {
    matched <- Filter(function(y) y$K == x$K && column %in% names(y$d), cases)
    values <- unlist(lapply(matched, function(y) c(y$d[[column]], theoretical(nrow(y$d), df))))
    if (is.null(df)) {
      limit <- ceiling(max(abs(values)) * 1.04 * 2) / 2
      c(-limit, limit)
    } else c(0, ceiling(max(values) * 1.04))
  }

  draw_qq <- function(x, column, statistic, df=NULL, node=FALSE) {
    plot_style()
    limits <- qq_limits(x, column, df)
    plot.new()
    plot.window(limits, limits, asp=1)
    abline(0,1,lty=2,lwd=0.85,col="#777777")
    points(theoretical(nrow(x$d),df),sort(x$d[[column]]),
           pch=1,cex=0.48,lwd=0.65,col="#285D89")
    draw_axes()
    add_labels(x,
      paste0(if (is.null(df)) "Normal Q-Q plot of $" else "Chi-square Q-Q plot of $", statistic, "$"),
      if (is.null(df)) "Standard normal quantiles" else paste0("$\\chi^2_{",df,"}$ quantiles"),
      "Empirical quantiles",node=node)
  }
  hist_pair <- function(x, j) {
    raw <- x$raw_error[, j]
    corrected <- x$debiased_error[, j]
    limits <- range(c(raw, corrected, 0))
    padding <- max(diff(limits) * 1e-8, 1e-12)
    breaks <- seq(limits[1] - padding, limits[2] + padding, length.out = n_bins + 1L)
    h_raw <- hist(raw, breaks = breaks, plot = FALSE, include.lowest = TRUE)
    h_corrected <- hist(corrected, breaks = breaks, plot = FALSE, include.lowest = TRUE)
    list(raw = h_raw, corrected = h_corrected, truth = x$true_entries[j],
         a = x$coordinates$a[j], b = x$coordinates$b[j])
  }

  draw_hist <- function(x, pair, errors=FALSE) {
    plot_style()
    shift <- if (errors) 0 else pair$truth
    breaks <- pair$raw$breaks + shift
    padding <- diff(range(breaks))*0.025
    peak <- max(pair$raw$density,pair$corrected$density)
    plot.new()
    plot.window(range(breaks)+c(-padding,padding),c(0,peak*1.25))
    rect(head(breaks,-1L),0,tail(breaks,-1L),pair$raw$density,
         col=raw_fill,border="#6D3831",lwd=0.55)
    rect(head(breaks,-1L),0,tail(breaks,-1L),pair$corrected$density,
         col=debiased_fill,border="#2C4269",lwd=0.55)
    abline(v=shift,col="#555555",lty=2,lwd=0.8)
    draw_axes()
    ij <- if (max(pair$a,pair$b)<10L) paste0(pair$a,pair$b) else paste(pair$a,pair$b,sep=",")
    estimate_names <- paste0(c("\\widehat{\\Lambda}_{","\\widetilde{\\Lambda}_{"),ij,"}")
    legend_names <- if (errors) paste0(estimate_names,"-\\Lambda^{\\star}_{",ij,"}") else estimate_names
    # Reserve legend space with plotmath; render labels in LaTeX.
    geometry_labels <- if (errors) as.expression(list(
      bquote(hat(Lambda)[.(ij)]-Lambda[.(ij)]^"*"),
      bquote(tilde(Lambda)[.(ij)]-Lambda[.(ij)]^"*"))) else as.expression(list(
      bquote(hat(Lambda)[.(ij)]),bquote(tilde(Lambda)[.(ij)])))
    geometry <- legend("topright",legend=geometry_labels,
      fill=c(raw_fill,debiased_fill),border=c("#6D3831","#2C4269"),
      bty="n",cex=1.12,inset=c(0.015,0.01),y.intersp=1.40,x.intersp=0.65,plot=FALSE)
    char_width <- strwidth("M",units="user",cex=1.12)
    box_width <- 0.8*char_width
    for (i in seq_len(2L)) {
      right <- geometry$text$x[i]-0.65*char_width
      left <- right-box_width
      h <- strheight("M",units="user",cex=1.12)*0.8
      rect(left,geometry$text$y[i]-h/2,right,geometry$text$y[i]+h/2,
           col=c(raw_fill,debiased_fill)[i],border=c("#6D3831","#2C4269")[i],lwd=0.7)
      label_at(xbp(geometry$text$x[i]),ybp(geometry$text$y[i]),
               paste0("$",legend_names[i],"$"),size=12.32,align="l")
    }
    heading <- if (errors) "Sampling distributions of estimation errors" else
      paste0("Sampling distributions of $",estimate_names[1],"$ and $",estimate_names[2],"$")
    add_labels(x,heading,if (errors) "Estimation error" else "Estimate","Density")
  }
  write_tex <- function(name) {
    commands <- vapply(labels, function(a) {
      content <- sprintf("{\\fontsize{%.3fbp}{%.3fbp}\\selectfont %s}",a$size,a$size*1.2,a$text)
      if (a$angle != 0) content <- sprintf("\\rotatebox{%d}{%s}",a$angle,content)
      sprintf("\\put(%.5f,%.5f){\\makebox(0,0)[%s]{%s}}",a$x,a$y,a$align,content)
    }, character(1))
    tex <- c(
      "\\documentclass{article}",
      "\\usepackage[paperwidth=5.5in,paperheight=5.5in,margin=0pt]{geometry}",
      "\\usepackage[T1]{fontenc}",
      "\\usepackage{lmodern,amsmath,amssymb,graphicx,xcolor}",
      "\\definecolor{ink}{HTML}{222222}",
      "\\pagestyle{empty}",
      "\\setlength{\\parindent}{0pt}",
      "\\setlength{\\topskip}{0pt}",
      "\\setlength{\\unitlength}{1bp}",
      "\\begin{document}",
      "\\color{ink}",
      "\\noindent\\begin{picture}(396,396)",
      paste0("\\put(0,0){\\makebox(0,0)[lb]{\\includegraphics[width=396bp,height=396bp]{",name,"_background.pdf}}}"),
      commands,
      "\\end{picture}",
      "\\end{document}")
    writeLines(tex,file.path(build_dir,paste0(name,".tex")))
  }
  compile_tex <- function(name, pdf_path) {
    engine <- Sys.which("pdflatex")
    if (!nzchar(engine)) stop("Install TinyTeX/TeX Live with pdflatex before rendering.")
    prior <- setwd(build_dir); on.exit(setwd(prior),add=TRUE)
    log <- system2(engine, c("-interaction=nonstopmode","-halt-on-error",
      "-file-line-error",shQuote(paste0(name,".tex"))),stdout=TRUE,stderr=TRUE)
    status <- attr(log,"status")
    if (!is.null(status) && status != 0) stop(paste(tail(log,25),collapse="\n"))
    ok <- file.copy(paste0(name,".pdf"),pdf_path,overwrite=TRUE)
    if (!ok) stop("Could not save final PDF.")
  }
  save_figure <- function(x,kind,component,a=NA_integer_,b=NA_integer_,draw) {
    folder <- if (kind %in% c("hist_estimates","hist_errors")) kind else "qq"
    tag <- if (is.na(component)) "joint" else if (!is.na(a)) {
      paste0("lambda_",a,"_",b)
    } else paste0("component",component)
    stem <- paste0(x$model,"_",x$scenario,"_n",x$cfg$n,"_",kind,"_",tag)
    pdf_path <- file.path(output_dir,folder,paste0(stem,".pdf"))
    png_path <- file.path(output_dir,"previews",paste0(stem,".png"))
    labels <<- list()
    cairo_pdf(file.path(build_dir,paste0(stem,"_background.pdf")),
      width=5.5,height=5.5,pointsize=11,family="serif")
    tryCatch(draw(),finally=dev.off())
    write_tex(stem)
    compile_tex(stem,pdf_path)
    if (nzchar(preview_engine)) {
      status <- system2(preview_engine,c("-singlefile","-r","150","-png",
        shQuote(pdf_path),shQuote(sub("[.]png$","",png_path))),stdout=FALSE,stderr=FALSE)
      if (status != 0L || !file.exists(png_path)) stop("Could not render PDF preview: ",pdf_path)
    }
    figure_records[[length(figure_records)+1L]] <<- data.frame(
      model=x$model,scenario=x$scenario,n=x$cfg$n,K=x$K,kind=kind,
      component=component,a=a,b=b,replications=nrow(x$d),
      pdf=file.path(folder,paste0(stem,".pdf")),
      preview=if (file.exists(png_path)) file.path("previews",paste0(stem,".png")) else NA_character_)
    cat(sprintf("[%03d] %s\n",length(figure_records),stem))
    flush.console()
  }

  if ("qq" %in% figure_types) {

  for (x in cases) {
    for (j in seq_len(x$K)) save_figure(x,"qq_theorem1",j,draw=function() {
      draw_qq(x,paste0("stat_U_",j),paste0("Z_{U,",j,"}"),node=TRUE)
    })
    for (j in seq_len(x$P)) save_figure(x,"qq_theorem2_debiased",j,draw=function() {
      draw_qq(x,paste0("debiased_lambda_",j),paste0("Z_{\\Lambda,",j,"}"))
    })
    save_figure(x,"qq_theorem1",NA_integer_,draw=function() {
      draw_qq(x,"chi_sq_U","\\|Z_U\\|^2",df=x$K,node=TRUE)
    })
    save_figure(x,"qq_theorem2_debiased",NA_integer_,draw=function() {
      draw_qq(x,"chi_sq_debiased_lambda","\\|Z_{\\Lambda}\\|^2",df=x$P)
    })
  }

  }

  if (any(c("hist_estimates", "hist_errors") %in% figure_types)) for (x in cases) for (j in seq_len(x$P)) {
    pair <- hist_pair(x, j)
    if ("hist_estimates" %in% figure_types) save_figure(x, "hist_estimates", j, pair$a, pair$b,
                draw = function() draw_hist(x, pair, errors = FALSE))
    if ("hist_errors" %in% figure_types) save_figure(x, "hist_errors", j, pair$a, pair$b,
                draw = function() draw_hist(x, pair, errors = TRUE))
    bin_records[[length(bin_records) + 1L]] <- data.frame(
      model = x$model, scenario = x$scenario, n = x$cfg$n, component = j,
      a = pair$a, b = pair$b, truth = pair$truth, bin = seq_len(n_bins),
      error_left = head(pair$raw$breaks, -1L), error_right = tail(pair$raw$breaks, -1L),
      estimate_left = head(pair$raw$breaks, -1L) + pair$truth,
      estimate_right = tail(pair$raw$breaks, -1L) + pair$truth,
      raw_count = pair$raw$counts, debiased_count = pair$corrected$counts,
      raw_density = pair$raw$density, debiased_density = pair$corrected$density)
    for (type in c("estimate", "error")) for (method in c("raw", "debiased")) {
      values <- x[[paste(method, type, sep = "_")]][, j]
      summary_records[[length(summary_records) + 1L]] <- data.frame(
        model = x$model, scenario = x$scenario, n = x$cfg$n, component = j,
        a = pair$a, b = pair$b, truth = pair$truth, quantity = type, method = method,
        replications = length(values), mean = mean(values), variance = var(values),
        min = min(values), max = max(values))
    }
  }
  figures <- do.call(rbind, figure_records)
  write.csv(figures, file.path(output_dir, "figure_manifest.csv"), row.names = FALSE)
  if (length(bin_records)) write.csv(do.call(rbind, bin_records), file.path(output_dir, "records", "histogram_bins.csv"), row.names = FALSE)
  if (length(summary_records)) write.csv(do.call(rbind, summary_records), file.path(output_dir, "records", "lambda_summary.csv"), row.names = FALSE)
  saveRDS(list(params = params, output_dir = output_dir, session = sessionInfo()),
          file.path(output_dir, "records", "plot_settings.rds"))
  expected <- sum(vapply(cases, function(x) {
    as.integer(("qq" %in% figure_types) * (x$K + x$P + 2L) +
      sum(c("hist_estimates", "hist_errors") %in% figure_types) * x$P)
  }, integer(1)))
  stopifnot(nrow(figures) == expected)
  stopifnot(!anyDuplicated(figures$pdf), all(file.exists(file.path(output_dir, figures$pdf))))
  cat("Output directory:", output_dir, "\n")

  input_hashes_unchanged <- vapply(seq_len(nrow(input_records)),function(i) {
    identical(unname(tools::md5sum(input_records$csv[i])),input_records$csv_md5[i]) &&
    identical(unname(tools::md5sum(input_records$config[i])),input_records$config_md5[i])
  },logical(1))
  stopifnot(all(input_hashes_unchanged))
  # Save a gallery with relative links.
  html <- c('<!doctype html><html lang="en"><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1">',
    paste0('<title>Individual figures at n=', n, '</title>'),
    '<style>body{max-width:1400px;margin:40px auto;padding:0 24px;color:#222;font-family:Georgia,serif}h1{font-weight:400}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:24px}figure{margin:0;border:1px solid #ddd;padding:12px}img{width:100%;display:block}figcaption{font:14px/1.5 system-ui,sans-serif}a{color:#285d89}</style>',
    paste0('<h1>Individual figures at n = ', n, '</h1>'),
    '<p>Figures use the saved cases listed in records/inputs.csv. Theorem 2 QQ plots show debiased statistics only. Text and mathematics are typeset in Latin Modern.</p>')
  for (kind in unique(figures$kind)) {
    html <- c(html,paste0('<h2>',kind,'</h2><div class="grid">'))
    for (i in which(figures$kind==kind)) {
      f <- figures[i,]
      label <- paste(f$model,f$scenario,
        if (is.na(f$component)) 'joint' else if (!is.na(f$a)) paste0('entry (',f$a,',',f$b,')') else paste('component',f$component),sep=' / ')
      html <- c(html,paste0('<figure><a href="',f$pdf,'">',
        if (!is.na(f$preview)) paste0('<img loading="lazy" src="',f$preview,'" alt="',label,'">') else 'Open PDF',
        '</a><figcaption>',label,' &mdash; <a href="',f$pdf,'">PDF</a></figcaption></figure>'))
    }
    html <- c(html,'</div>')
  }
  writeLines(c(html,'</html>'),file.path(output_dir,'figure_gallery.html'))

  return(invisible(list(output_dir = output_dir, figures = figures, inputs = input_records)))
}
