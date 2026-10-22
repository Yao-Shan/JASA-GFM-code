# Plot saved means and sample variances; sourcing only defines the function.
gfb_plot_summary <- function(input, output) {
  if (!file.exists(input)) stop("Summary CSV does not exist.")
  if (file.exists(output)) stop("Output already exists; choose a new PDF name.")
  d <- read.csv(input, stringsAsFactors = FALSE)
  needed <- c("scenario", "model", "theorem", "estimator", "coordinate", "mean", "sample_variance", "valid_n")
  if (!all(needed %in% names(d))) stop("Unsupported summary schema; missing: ", paste(setdiff(needed, names(d)), collapse = ", "))
  if (any(!is.finite(d$valid_n) | d$valid_n < 0)) stop("Invalid valid_n.")
  groups <- unique(d[, c("scenario", "theorem", "estimator"), drop = FALSE])
  grDevices::pdf(output, width = 10, height = 6, family = "Times", onefile = TRUE, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (i in seq_len(nrow(groups))) {
    g <- groups[i, ]
    x <- d[d$scenario == g$scenario & d$theorem == g$theorem & d$estimator == g$estimator, , drop = FALSE]
    x <- x[order(x$model, x$coordinate), , drop = FALSE]
    labels <- paste(x$model, x$coordinate, sep = " / ")
    for (stat in c("mean", "sample_variance")) {
      values <- x[[stat]]; ref <- if (stat == "mean") 0 else 1
      valid <- is.finite(values) & x$valid_n > 0
      par(mar = c(9, 4, 3, 1), family = "serif", mgp = c(2.5, .7, 0))
      heading <- paste(g$theorem, g$estimator, g$scenario, sep = ", ")
      if (!any(valid)) {plot.new(); title(main = heading); text(.5, .5, "No finite statistics"); next}
      limits <- range(c(values[valid], ref)); if (diff(limits) == 0) limits <- limits + c(-.05, .05)
      plot(which(valid), values[valid], pch = 1, col = "#285D89", xaxt = "n", xlab = "",
           xlim = c(.5, nrow(x) + .5), ylab = if(stat == "mean") "Mean" else "Sample variance",
           ylim = limits, main = heading, font.main = 1, bty = "l")
      axis(1, at = seq_along(values), labels = labels, las = 2, cex.axis = .7)
      abline(h = ref, lty = 2, col = "gray50")
    }
  }
  invisible(output)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) != 2L) stop("Usage: Rscript plots/plot_summary.R data/summary.csv output.pdf")
  gfb_plot_summary(args[1L], args[2L])
  cat("Saved", args[2L], "from", args[1L], "; no fitting.\n")
}
