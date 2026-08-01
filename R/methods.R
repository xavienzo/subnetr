#' @export
print.subnet_partition <- function(x, ...) {
  cat("Greedy-peeling partition\n")
  cat(sprintf("  nodes           : %d\n", x$n_nodes))
  cat(sprintf("  objective       : %s (lambda = %.2f)\n", x$objective, x$lambda))
  cat(sprintf("  threshold       : %.4g  (%.1f%% of edges retained)\n",
              x$threshold, 100 * x$pi0))
  cat(sprintf("  subnetworks     : %d\n", x$n_clusters))
  if (x$n_clusters > 0L) {
    i <- seq_len(x$n_clusters)
    df <- data.frame(subnet = i, size = x$sizes[i], edges = x$n_edge[i],
                     density = round(x$density[i], 3))
    cat("\n")
    print(df, row.names = FALSE)
  }
  bg <- length(x$sizes)
  cat(sprintf("\n  background      : %d nodes\n", x$sizes[bg]))
  invisible(x)
}

#' @export
print.subnet <- function(x, ...) {
  cat("Predictor-associated subnetworks\n\n")
  cat(sprintf("  nodes       : %d   edges: %d\n", x$partition$n_nodes,
              x$partition$n_nodes * (x$partition$n_nodes - 1L) / 2L))
  cat(sprintf("  objective   : %s (lambda = %.2f)\n", x$objective, x$lambda))
  cat(sprintf("  threshold   : %.4g   (%.1f%% of edges retained)\n",
              x$threshold, 100 * x$pi0))
  cat(sprintf("  permutations: %d   alpha: %.3g\n", x$n_perm, x$alpha))
  if (!is.null(x$tuning)) {
    cat(sprintf("  tuning      : %s criterion\n", x$tuning$criterion))
  }
  cat("\n")

  if (length(x$size) == 0L) {
    cat("No subnetwork was extracted.\n")
    return(invisible(x))
  }

  df <- as.data.frame(x)
  df$density <- round(df$density, 3)
  df$log10_stat <- round(df$statistic / log(10), 1)
  df$statistic <- NULL
  df$p_value <- format.pval(df$p_value, digits = 2, eps = 1 / (x$n_perm + 1))
  print(df, row.names = FALSE)

  ns <- sum(x$significant)
  cat(sprintf("\n%d of %d subnetwork%s significant at FWER %.3g.\n",
              ns, length(x$size), if (length(x$size) == 1L) "" else "s",
              x$alpha))
  cat("P-values are max-statistic permutation p-values and are already\n",
      "family-wise-error corrected; do not adjust them again.\n", sep = "")
  invisible(x)
}

#' Tabulate detected subnetworks
#'
#' @param x A `subnet` object.
#' @param row.names,optional Ignored, present for method consistency.
#' @param ... Ignored.
#'
#' @return A data frame with one row per extracted subnetwork: its index,
#'   node count, supra-threshold edge count, within-block density, log binomial
#'   statistic, permutation p-value and significance flag.
#'
#' @examples
#' sim <- simulate_fc(n = 100, n_nodes = 60, cluster_size = 12, f2 = 0.2,
#'                    seed = 5)
#' as.data.frame(subnet(sim$W, n_perm = 99, seed = 1))
#'
#' @export
as.data.frame.subnet <- function(x, row.names = NULL, optional = FALSE, ...) {
  data.frame(subnet = seq_along(x$size),
             size = x$size,
             edges = x$n_edge,
             density = x$density,
             statistic = x$statistic,
             p_value = x$p_value,
             significant = x$significant)
}

#' Subnetwork membership for every node
#'
#' @param x A `subnet` object.
#' @param significant_only Label only the subnetworks that reached
#'   significance, leaving the rest as background.
#'
#' @return An integer vector of length `n_nodes`. Entry `i` is the index of the
#'   subnetwork containing node `i`, or `0` for background nodes.
#'
#' @examples
#' sim <- simulate_fc(n = 100, n_nodes = 60, cluster_size = 12, f2 = 0.2,
#'                    seed = 5)
#' table(membership(subnet(sim$W, n_perm = 99, seed = 1)))
#'
#' @export
membership <- function(x, significant_only = TRUE) {
  stopifnot(inherits(x, "subnet"))
  out <- integer(x$partition$n_nodes)
  keep <- if (significant_only) which(x$significant) else seq_along(x$size)
  for (i in keep) out[x$subnetworks[[i]]] <- i
  out
}

#' @export
print.subnet_tuning <- function(x, ...) {
  cat("Threshold / lambda tuning\n")
  cat(sprintf("  criterion : %s%s\n", x$criterion,
              if (x$criterion == "calibrated")
                sprintf(" (%d permutations per grid point)", x$n_perm) else ""))
  cat(sprintf("  grid      : %d lambda x %d threshold\n",
              length(unique(x$grid$lambda)),
              length(unique(x$grid$threshold))))
  cat(sprintf("  selected  : lambda = %.2f, threshold = %.4g\n",
              x$lambda, x$threshold))
  cat("\nTop grid points:\n")
  score <- if (x$criterion == "calibrated") x$grid$z else x$grid$lr
  ord <- order(score, decreasing = TRUE)[seq_len(min(5L, nrow(x$grid)))]
  show <- x$grid[ord, c("lambda", "prob", "threshold", "n_clusters", "lr", "z")]
  show$threshold <- signif(show$threshold, 4)
  show$lr <- round(show$lr, 1)
  show$z <- round(show$z, 2)
  print(show, row.names = FALSE)
  invisible(x)
}

#' @export
print.subnet_power <- function(x, ...) {
  p <- x$params
  cat("Power analysis for subnetwork detection\n\n")
  cat(sprintf("  nodes           : %d\n", p$n_nodes))
  cat(sprintf("  planted subnets : %s (sizes)\n",
              paste(p$cluster_size, collapse = ", ")))
  cat(sprintf("  effect size f2  : %s\n", paste(p$f2, collapse = ", ")))
  cat(sprintf("  rho_in / rho_out: %.2f / %.3f\n", p$rho_in, p$rho_out))
  cat(sprintf("  replicates      : %d   permutations: %d   alpha: %.3g\n",
              p$n_sim, p$n_perm, p$alpha))
  cat(sprintf("  Dice cutoff     : %.2f\n\n", p$dice_cut))

  d <- x$power
  show <- data.frame(
    n = d$n,
    power = sprintf("%.3f (%.3f)", d$power, d$se_power),
    recovery = sprintf("%.3f (%.3f)", d$power_recovery, d$se_recovery),
    dice = round(d$mean_dice, 3),
    n_sig = round(d$mean_n_sig, 2))
  print(show, row.names = FALSE)
  cat("\nPower and recovery show the estimate with its Monte Carlo standard\n",
      "error in parentheses.\n", sep = "")
  invisible(x)
}

#' @export
print.subnet_sim <- function(x, ...) {
  p <- x$params
  cat("Simulated connectome study\n")
  cat(sprintf("  subjects  : %d   nodes: %d   covariates: %d\n",
              p$n, p$n_nodes, p$n_cov))
  cat(sprintf("  planted   : %s nodes per subnetwork\n",
              paste(p$cluster_size, collapse = ", ")))
  cat(sprintf("  f2        : %s   rho_in: %.2f   rho_out: %.3f\n",
              paste(p$f2, collapse = ", "), p$rho_in, p$rho_out))
  cat(sprintf("  method    : %s   signal edges: %d of %d\n",
              p$method, sum(x$truth$edge_signal), length(x$truth$edge_signal)))
  invisible(x)
}

## Heatmap of a square matrix in conventional matrix orientation, i.e. row 1
## at the top.
image_matrix <- function(M, col, zlim, main = NULL, xlab = "node",
                         ylab = "node", ...) {
  n <- nrow(M)
  graphics::image(x = seq_len(n), y = seq_len(n), z = t(M[n:1, , drop = FALSE]),
                  col = col, zlim = zlim, axes = FALSE, main = main,
                  xlab = xlab, ylab = ylab, ...)
  # Rows are drawn top-down, so a row plotted at height p is node n - p + 1.
  labs <- pretty(seq_len(n))
  labs <- labs[labs >= 1 & labs <= n]
  graphics::axis(1, at = labs, labels = labs)
  graphics::axis(2, at = n - labs + 1, labels = labs, las = 1)
  graphics::box()
}

## Vertical colour scale annotated with the quantity being displayed. Drawn in
## its own plot region, so the caller is responsible for the layout.
color_bar <- function(zlim, col, label) {
  op <- graphics::par(mar = c(5.1, 0.6, 4.1, 3.6))
  on.exit(graphics::par(op), add = TRUE)
  z <- seq(zlim[1], zlim[2], length.out = length(col) + 1L)
  graphics::image(x = 1, y = z, z = matrix(z[-length(z)], nrow = 1),
                  col = col, zlim = zlim, axes = FALSE, xlab = "", ylab = "")
  graphics::axis(4, las = 1, cex.axis = 0.8, mgp = c(3, 0.5, 0), tcl = -0.25)
  graphics::box()
  if (!is.null(label)) {
    graphics::mtext(label, side = 3, line = 0.6, cex = 0.85, adj = 0.5)
  }
}

#' Plot the connectivity matrix
#'
#' Draws the edge-weight matrix, optionally with nodes reordered so that the
#' extracted subnetworks sit in consecutive blocks along the diagonal, and
#' outlines the significant blocks.
#'
#' `what = "both"` puts the matrix as measured next to the reordered version,
#' sharing one colour scale. That pairing is the honest way to present a
#' result: a subnetwork is invisible in the original node ordering and obvious
#' after reordering, and showing only the second panel can make an arbitrary
#' permutation look like a discovery. Comparing the two is also how you judge
#' whether a detection is a compact block or a diffuse smear.
#'
#' The colour scale is labelled with the quantity being displayed. When `W`
#' came from [edge_stats()] or [simulate_fc()] that label is carried along with
#' the matrix, so the figure states whether it is showing \eqn{-\log_{10}} p,
#' a t statistic, or something else. Supply `weight_label` for matrices built
#' by other means.
#'
#' @param x A `subnet` object.
#' @param what `"reordered"` (default), `"observed"` for the original node
#'   order, or `"both"` for the two side by side.
#' @param significant_only Outline only the significant subnetworks.
#' @param col Colour palette.
#' @param main Panel title, or a vector of two when `what = "both"`.
#' @param weight_label Label for the colour scale, as a string or an
#'   [expression]. Defaults to the `"statistic"` attribute of `W`, and to
#'   `"edge weight"` when the matrix carries none.
#' @param legend Draw the colour scale.
#' @param ... Passed to [graphics::image()].
#'
#' @return `x`, invisibly. Called for the plot.
#'
#' @examples
#' sim <- simulate_fc(n = 120, n_nodes = 60, cluster_size = 12, f2 = 0.2,
#'                    seed = 11)
#' fit <- subnet(sim$W, n_perm = 99, seed = 1)
#'
#' plot(fit)                  # reordered only
#' plot(fit, what = "both")   # as measured, next to reordered
#'
#' # A matrix of t statistics from a real study labels itself accordingly
#' plot(fit, weight_label = "t statistic")
#'
#' @export
plot.subnet <- function(x, what = c("reordered", "observed", "both"),
                        significant_only = TRUE,
                        col = grDevices::hcl.colors(64, "Inferno"),
                        main = NULL, weight_label = NULL, legend = TRUE, ...) {
  what <- match.arg(what)
  W <- x$W
  # A single scale across both panels: they show the same numbers, so a shared
  # range is what makes them comparable.
  zlim <- c(0, max(W))
  if (is.null(weight_label)) {
    weight_label <- attr(W, "statistic") %||% "edge weight"
  }

  panels <- if (what == "both") c("observed", "reordered") else what
  titles <- if (is.null(main)) {
    c(observed = "Original node order",
      reordered = "Reordered by subnetwork")[panels]
  } else {
    rep_len(main, length(panels))
  }

  if (legend || length(panels) > 1L) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit({graphics::layout(1); graphics::par(op)}, add = TRUE)
    widths <- c(rep(5, length(panels)), if (legend) 1.35)
    graphics::layout(matrix(seq_along(widths), nrow = 1), widths = widths)
    if (length(panels) > 1L) graphics::par(mar = c(4.6, 4.0, 3.2, 1.0))
  }

  for (i in seq_along(panels)) {
    if (panels[i] == "observed") {
      image_matrix(W, col, zlim, main = titles[i], ...)
      next
    }
    ord <- x$partition$order
    n <- length(ord)
    image_matrix(W[ord, ord, drop = FALSE], col, zlim, main = titles[i], ...)

    keep <- if (significant_only) which(x$significant) else seq_along(x$size)
    if (length(keep)) {
      ends <- cumsum(x$partition$sizes)
      starts <- ends - x$partition$sizes + 1L
      for (k in keep) {
        # image_matrix flips the vertical axis, so the y range is mirrored.
        graphics::rect(starts[k] - 0.5, n - ends[k] + 0.5,
                       ends[k] + 0.5, n - starts[k] + 1.5,
                       border = "#39FF14", lwd = 2.5)
      }
    }
  }

  if (legend) color_bar(zlim, col, weight_label)
  invisible(x)
}

#' Plot a power curve
#'
#' @param x A `subnet_power` object.
#' @param which Curves to draw: `"both"` (default), `"any"` or `"recovery"`.
#' @param ci Draw pointwise Monte Carlo 95% intervals.
#' @param target Optional horizontal reference line, e.g. `0.8`.
#' @param ... Passed to [graphics::plot()].
#'
#' @return `x`, invisibly. Called for the plot.
#'
#' @examples
#' \donttest{
#' pw <- power_curve(n = c(40, 80, 120), n_nodes = 60, cluster_size = 12,
#'                   f2 = 0.2, n_sim = 20, n_perm = 99, seed = 1)
#' plot(pw, target = 0.8)
#' }
#'
#' @export
plot.subnet_power <- function(x, which = c("both", "any", "recovery"),
                              ci = TRUE, target = 0.8, ...) {
  which <- match.arg(which)
  d <- x$power

  graphics::plot(range(d$n), c(0, 1), type = "n",
                 xlab = "sample size", ylab = "power",
                 main = "Subnetwork detection power", ...)
  graphics::grid(col = "grey90", lty = 1)
  if (!is.null(target)) {
    graphics::abline(h = target, lty = 2, col = "grey40")
  }

  draw <- function(y, se, colr, pch, lty, lwd) {
    if (ci) {
      lo <- pmax(0, y - 1.96 * se); hi <- pmin(1, y + 1.96 * se)
      # A proportion of exactly 0 or 1 has zero standard error; drawing a
      # zero-length arrow warns and renders nothing useful.
      k <- which(hi - lo > .Machine$double.eps)
      if (length(k)) {
        graphics::arrows(d$n[k], lo[k], d$n[k], hi[k], angle = 90, code = 3,
                         length = 0.03, col = colr)
      }
    }
    graphics::lines(d$n, y, col = colr, lwd = lwd, lty = lty)
    graphics::points(d$n, y, col = colr, pch = pch, bg = "white")
  }

  # Recovery frequently equals detection once the effect is easy to find, and
  # is drawn second. A dashed overlay on a thicker solid line keeps both
  # readable where they coincide instead of one hiding the other.
  leg <- character(0); cols <- character(0)
  pchs <- numeric(0); ltys <- numeric(0)
  if (which %in% c("both", "any")) {
    draw(d$power, d$se_power, "#1f77b4", 21, 1, 3)
    leg <- c(leg, "any detection"); cols <- c(cols, "#1f77b4")
    pchs <- c(pchs, 21); ltys <- c(ltys, 1)
  }
  if (which %in% c("both", "recovery")) {
    draw(d$power_recovery, d$se_recovery, "#d62728", 22, 2, 1.8)
    leg <- c(leg, sprintf("recovery (Dice >= %.2f)", x$params$dice_cut))
    cols <- c(cols, "#d62728"); pchs <- c(pchs, 22); ltys <- c(ltys, 2)
  }
  graphics::legend("bottomright", legend = leg, col = cols, pch = pchs,
                   lty = ltys, lwd = c(3, 1.8)[seq_along(leg)], bty = "n")
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
