#' Simulate a connectome study with planted predictor-associated subnetworks
#'
#' Generates a weighted adjacency matrix of \eqn{-\log_{10}} p-values in which
#' one or more subnetworks carry an association with a predictor of interest,
#' together with the ground truth needed to score a detection method.
#'
#' @details
#' Nodes 1 through `sum(cluster_size)` are assigned to consecutive planted
#' subnetworks. A proportion `rho_in` of the edges inside each subnetwork
#' actually carry the effect (the rest are within-subnetwork nulls, reflecting
#' that a real subnetwork is not uniformly affected), and a proportion
#' `rho_out` of the edges outside every subnetwork carry it too (scattered
#' associations that are not part of any subnetwork). Node labels are then
#' shuffled so the planted structure is not sitting in the corner of the matrix.
#'
#' @section Simulation method:
#' `method = "fast"` (default) draws the edge t-statistics directly from their
#' exact sampling distribution rather than materializing subject-level data.
#' Under the linear model with fixed design, an edge with Cohen's
#' \eqn{f^2} effect size has
#' \deqn{t \sim t_{df}(\mathrm{ncp}), \quad \mathrm{ncp} = \sqrt{f^2 S},
#'   \quad S \sim \chi^2_{n - 1 - q},}
#' where \eqn{df = n - 2 - q} and \eqn{q} is the number of nuisance covariates.
#' `S` is drawn once per data set, which is what correctly couples all edges
#' through the shared predictor. This is exact, not an approximation, and it
#' avoids allocating an `n` by `n_edge` matrix -- the reason power analysis
#' over thousands of replicates is practical.
#'
#' `method = "data"` builds the subject-level edge matrix explicitly and runs
#' [edge_stats()] on it. Use it when you need the raw data, or when
#' `edge_corr > 0` to study how the test behaves under cross-edge dependence,
#' which the fast path deliberately excludes.
#'
#' @param n Number of subjects.
#' @param n_nodes Number of nodes in the connectome.
#' @param cluster_size Integer vector of planted subnetwork sizes.
#' @param f2 Cohen's \eqn{f^2} effect size for the affected edges, recycled
#'   across subnetworks so different subnetworks can carry different effects.
#' @param rho_in Proportion of within-subnetwork edges that carry the effect.
#' @param rho_out Proportion of outside edges that carry the effect.
#' @param f2_out Effect size for the scattered outside edges. Defaults to
#'   `f2[1]`.
#' @param n_cov Number of nuisance covariates adjusted for.
#' @param method `"fast"` or `"data"`; see the section on simulation method.
#' @param edge_corr Cross-edge noise correlation induced by a shared
#'   subject-level factor. `method = "data"` only.
#' @param shuffle_nodes Randomly relabel nodes so the planted subnetworks are
#'   scattered through the matrix.
#' @param seed Optional integer seed.
#'
#' @return A list of class `subnet_sim`:
#'   \item{W}{`n_nodes` by `n_nodes` matrix of \eqn{-\log_{10}} p-values.}
#'   \item{truth}{List with `nodes` (node indices of each planted subnetwork,
#'     in the returned labelling), `edge_signal` (logical over vectorized
#'     edges), and `labels` (block index per node, 0 for background).}
#'   \item{data}{Subject-level `fc`, `x` and `covariates`, or `NULL` under
#'     `method = "fast"`.}
#'   \item{params}{The settings used, including the residual `df`.}
#'
#' @seealso [subnet()], [power_curve()], [edge_stats()]
#'
#' @examples
#' sim <- simulate_fc(n = 100, n_nodes = 80, cluster_size = c(20, 10),
#'                    f2 = c(0.15, 0.25), seed = 3)
#' dim(sim$W)
#' lengths(sim$truth$nodes)
#'
#' @export
simulate_fc <- function(n, n_nodes = 100L, cluster_size = 20L, f2 = 0.05,
                        rho_in = 0.9, rho_out = 0.02, f2_out = NULL,
                        n_cov = 0L, method = c("fast", "data"),
                        edge_corr = 0, shuffle_nodes = TRUE, seed = NULL) {
  method <- match.arg(method)
  if (!is.null(seed)) set.seed(seed)

  n <- as.integer(n)
  n_nodes <- as.integer(n_nodes)
  n_cov <- as.integer(n_cov)
  cluster_size <- as.integer(cluster_size)

  if (any(cluster_size < 2L)) {
    stop("Every entry of `cluster_size` must be at least 2.", call. = FALSE)
  }
  if (sum(cluster_size) > n_nodes) {
    stop("`cluster_size` sums to ", sum(cluster_size),
         ", which exceeds `n_nodes` (", n_nodes, ").", call. = FALSE)
  }
  df <- n - 2L - n_cov
  if (df < 1L) {
    stop("`n` must exceed `n_cov` + 2 to leave residual degrees of freedom.",
         call. = FALSE)
  }
  f2 <- rep_len(f2, length(cluster_size))
  if (is.null(f2_out)) f2_out <- f2[1]

  M <- n_nodes * (n_nodes - 1L) / 2L

  # --- planted structure -----------------------------------------------------
  labels <- integer(n_nodes)
  ends <- cumsum(cluster_size)
  starts <- ends - cluster_size + 1L
  for (i in seq_along(cluster_size)) labels[starts[i]:ends[i]] <- i

  lab_mat <- matrix(0L, n_nodes, n_nodes)
  for (i in seq_along(cluster_size)) {
    lab_mat[starts[i]:ends[i], starts[i]:ends[i]] <- i
  }
  edge_block <- lab_mat[lower.tri(lab_mat)]

  # Which edges actually carry the effect: exact counts, sampled without
  # replacement, so replicate-to-replicate variance comes from the data rather
  # than from the contamination rate.
  signal <- logical(M)
  f2_edge <- numeric(M)
  for (i in seq_along(cluster_size)) {
    idx <- which(edge_block == i)
    take <- sample(idx, floor(rho_in * length(idx)))
    signal[take] <- TRUE
    f2_edge[take] <- f2[i]
  }
  out_idx <- which(edge_block == 0L)
  if (rho_out > 0 && length(out_idx)) {
    take <- sample(out_idx, floor(rho_out * length(out_idx)))
    signal[take] <- TRUE
    f2_edge[take] <- f2_out
  }

  # --- edge statistics -------------------------------------------------------
  dat <- NULL
  if (method == "fast") {
    S <- rchisq(1, df = n - 1L - n_cov)
    tstat <- numeric(M)
    tstat[!signal] <- rt(sum(!signal), df = df)
    for (val in unique(f2_edge[signal])) {
      j <- which(signal & f2_edge == val)
      tstat[j] <- rt(length(j), df = df, ncp = sqrt(val * S))
    }
    nlp <- -(log(2) + pt(-abs(tstat), df, log.p = TRUE)) / log(10)
    W <- unvech(nlp)
    attr(W, "statistic") <- expression(-log[10](p))
  } else {
    x <- rnorm(n)
    Z <- if (n_cov > 0L) matrix(rnorm(n * n_cov), n, n_cov) else NULL

    R2 <- f2_edge / (1 + f2_edge)
    E <- matrix(rnorm(n * M), n, M)
    if (edge_corr > 0) {
      shared <- rnorm(n)
      E <- sqrt(1 - edge_corr) * E + sqrt(edge_corr) * shared
    }
    Y <- E * rep(sqrt(1 - R2), each = n) + outer(x, sqrt(R2))
    if (!is.null(Z)) {
      G <- matrix(rnorm(n_cov * M, sd = 0.5), n_cov, M)
      Y <- Y + Z %*% G
    }
    W <- edge_stats(Y, x, covariates = Z, value = "nlog10p")
    dat <- list(fc = Y, x = x, covariates = Z)
  }

  # --- relabel nodes ---------------------------------------------------------
  perm <- if (shuffle_nodes) sample.int(n_nodes) else seq_len(n_nodes)
  lab <- attr(W, "statistic")
  W <- W[perm, perm, drop = FALSE]   # subsetting drops attributes
  attr(W, "statistic") <- lab
  labels <- labels[perm]
  # `perm` maps new position -> old node, so the new index of old node j is
  # order(perm)[j].
  inv <- order(perm)
  nodes <- lapply(seq_along(cluster_size), function(i) sort(inv[starts[i]:ends[i]]))

  # Re-express the signal indicator in the new labelling.
  sig_mat <- matrix(0L, n_nodes, n_nodes)
  sig_mat[lower.tri(sig_mat)] <- as.integer(signal)
  sig_mat <- sig_mat + t(sig_mat)
  edge_signal <- as.logical(vech(sig_mat[perm, perm, drop = FALSE]))

  structure(
    list(W = W,
         truth = list(nodes = nodes, edge_signal = edge_signal,
                      labels = labels),
         data = dat,
         params = list(n = n, n_nodes = n_nodes, cluster_size = cluster_size,
                       f2 = f2, f2_out = f2_out, rho_in = rho_in,
                       rho_out = rho_out, n_cov = n_cov, df = df,
                       method = method, edge_corr = edge_corr)),
    class = "subnet_sim")
}
