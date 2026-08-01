#' Power analysis for subnetwork detection
#'
#' Estimates, by simulation, the probability that a study of a given size
#' detects a predictor-associated subnetwork of a given size and effect size.
#'
#' @details
#' For each sample size and each replicate the function simulates a connectome
#' study with [simulate_fc()], runs the full [subnet()] pipeline including its
#' permutation test, and records what was found. Two quantities are reported
#' because they answer different questions:
#'
#' \describe{
#'   \item{`power`}{the probability of declaring *any* subnetwork significant.
#'     This is the conventional notion of power, but on its own it can be
#'     satisfied by a subnetwork that has little to do with the planted one.}
#'   \item{`power_recovery`}{the probability that every planted subnetwork is
#'     matched by a significant one with Sorensen-Dice overlap at least
#'     `dice_cut`. This is the quantity you want when the scientific claim is
#'     about *which* nodes are involved, and it is always the more demanding of
#'     the two.}
#' }
#'
#' The gap between them is informative: when `power` is high but
#' `power_recovery` is not, the study is large enough to see that something is
#' there but not to localize it.
#'
#' Replicates are independent, so the Monte Carlo standard error of each
#' proportion is reported and shrinks as `1 / sqrt(n_sim)`. Roughly 200
#' replicates give a standard error near 0.035 at 50% power.
#'
#' @section The screening threshold:
#' Nothing else in the analysis moves power as much as `threshold_prob`. A
#' subnetwork of `c` nodes spans `c(c-1)/2` edges; if screening retains fewer
#' edges than that across the whole graph, the subnetwork cannot survive as a
#' block no matter how large the sample. In an 80-node connectome with a
#' 20-node planted subnetwork at f2 = 0.06 and n = 120, recovery probability
#' runs
#'
#' \tabular{lr}{
#'   `threshold_prob` \tab recovery \cr
#'   0.99  \tab 0.30 \cr
#'   0.975 \tab 0.95 \cr
#'   0.95  \tab 1.00 \cr
#' }
#'
#' purely from where the threshold sits: the subnetwork spans 190 edges, and at
#' the 99th percentile only 32 edges survive screening in the entire graph. The
#' `NULL` default scales the threshold with the design to keep this from
#' silently dominating the answer, but a reported power is a statement about a
#' *particular* screening choice and should be quoted alongside it.
#'
#' If you intend to tune the threshold at analysis time rather than fix it in
#' advance, set `tune_each = TRUE` so the simulation runs the procedure you
#' will actually use. That is slower, and it can come out either higher or
#' lower than a fixed threshold -- permutations now pay the same search cost,
#' which costs power, but the threshold adapts to each data set, which gains
#' it. The reason to use it is fidelity, not conservatism.
#'
#' @param n Integer vector of sample sizes to evaluate.
#' @param n_nodes,cluster_size,f2,rho_in,rho_out,f2_out,n_cov Passed to
#'   [simulate_fc()], describing the connectome and the planted effect.
#' @param n_sim Replicates per sample size.
#' @param n_perm Permutations per replicate. `199` gives p-value resolution
#'   `0.005`, ample for `alpha = 0.05`.
#' @param alpha Family-wise significance level.
#' @param threshold_prob Quantile of the observed edge weights used as the
#'   screening threshold when `tune_each = FALSE`. `NULL` (the default) scales
#'   it with the design so that screening retains roughly twice as many edges
#'   as the planted subnetworks contain, capped at the top 10%. See the section
#'   on the screening threshold -- this is the single most consequential
#'   setting in a power analysis.
#' @param lambda Objective parameter used when `tune_each = FALSE`.
#' @param objective,min_size,max_clusters Passed to [subnet()].
#' @param tune_each Retune the threshold and `lambda` in every replicate with
#'   [tune_subnet()], and pay the corresponding selection cost in the
#'   permutation null. This simulates the procedure a real analysis runs, which
#'   is the reason to use it; it is considerably slower, and the power it
#'   reports can land either side of a fixed threshold.
#' @param tune A named list of extra arguments for [tune_subnet()], used only
#'   when `tune_each = TRUE`. Shrinking its grid or its `n_perm` is the main
#'   way to make `tune_each` affordable.
#' @param dice_cut Minimum Sorensen-Dice overlap for a planted subnetwork to
#'   count as recovered.
#' @param n_cores Cores used to run replicates in parallel.
#' @param seed Optional integer seed. Results are reproducible and do not
#'   depend on `n_cores`.
#' @param progress Print a line per sample size.
#'
#' @return An object of class `subnet_power`:
#'   \item{power}{Data frame with one row per sample size: `n`, `power`,
#'     `power_recovery`, their Monte Carlo standard errors `se_power` and
#'     `se_recovery`, `mean_dice`, `mean_n_sig` and `false_subnet_rate`, the
#'     average proportion of significant subnetworks that match no planted
#'     one.}
#'   \item{by_cluster}{Data frame with one row per sample size and planted
#'     subnetwork, giving that subnetwork's own recovery probability and mean
#'     overlap.}
#'   \item{replicates}{Per-replicate results, for diagnostics.}
#'   \item{params}{The settings used.}
#'
#' @seealso [required_n()], [plot.subnet_power()], [simulate_fc()], [subnet()]
#'
#' @examples
#' \donttest{
#' pw <- power_curve(n = c(50, 100), n_nodes = 60, cluster_size = 12,
#'                   f2 = 0.2, n_sim = 20, n_perm = 99, seed = 1)
#' pw
#' }
#'
#' @export
power_curve <- function(n, n_nodes = 100L, cluster_size = 20L, f2 = 0.05,
                        rho_in = 0.9, rho_out = 0.02, f2_out = NULL,
                        n_cov = 0L, n_sim = 200L, n_perm = 199L, alpha = 0.05,
                        threshold_prob = NULL, lambda = 0.6,
                        objective = c("generalized", "density"),
                        min_size = 3L, max_clusters = 25L,
                        tune_each = FALSE, tune = list(), dice_cut = 0.5,
                        n_cores = 1L, seed = NULL, progress = TRUE) {
  objective <- match.arg(objective)
  n <- sort(unique(as.integer(n)))
  n_sim <- as.integer(n_sim)
  n_perm <- as.integer(n_perm)
  cluster_size <- as.integer(cluster_size)
  n_cores <- max(1L, as.integer(n_cores))
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  n_clust <- length(cluster_size)

  # Screening must retain enough edges for the target subnetwork to survive as
  # a block. A threshold keeping fewer edges than the subnetwork contains
  # cannot possibly recover it, so the default scales with the design rather
  # than sitting at a fixed quantile.
  if (is.null(threshold_prob)) {
    m_signal <- sum(cluster_size * (cluster_size - 1) / 2)
    m_total <- n_nodes * (n_nodes - 1) / 2
    threshold_prob <- 1 - min(0.10, max(0.005, 2 * m_signal / m_total))
  }

  one_rep <- function(nn, rep_id, rep_seed) {
    set.seed(rep_seed)
    sim <- simulate_fc(n = nn, n_nodes = n_nodes, cluster_size = cluster_size,
                       f2 = f2, rho_in = rho_in, rho_out = rho_out,
                       f2_out = f2_out, n_cov = n_cov, method = "fast",
                       shuffle_nodes = TRUE)

    args <- list(W = sim$W, alpha = alpha, n_perm = n_perm,
                 objective = objective, min_size = min_size,
                 max_clusters = max_clusters, seed = rep_seed, n_cores = 1L)
    if (tune_each) {
      args$tune <- tune
    } else {
      args$threshold <- as.numeric(quantile(vech(sim$W), threshold_prob,
                                            names = FALSE))
      args$lambda <- lambda
    }
    fit <- do.call(subnet, args)

    sig <- which(fit$significant)
    # Best overlap between each planted subnetwork and any significant one.
    d <- vapply(sim$truth$nodes, function(tn) {
      if (!length(sig)) return(0)
      max(vapply(fit$subnetworks[sig], dice, numeric(1), b = tn))
    }, numeric(1))
    # Significant subnetworks matching no planted one.
    false_rate <- if (!length(sig)) NA_real_ else {
      mean(vapply(fit$subnetworks[sig], function(sn) {
        max(vapply(sim$truth$nodes, dice, numeric(1), a = sn)) < dice_cut
      }, logical(1)))
    }

    c(list(n = nn, rep = rep_id, detected = length(sig) > 0L,
           n_sig = length(sig),
           recovered_all = all(d >= dice_cut),
           mean_dice = mean(d),
           false_subnet_rate = false_rate),
      stats::setNames(as.list(d), paste0("dice", seq_len(n_clust))))
  }

  reps <- list()
  for (i in seq_along(n)) {
    if (progress) {
      message(sprintf("n = %d  (%d replicates)", n[i], n_sim))
      utils::flush.console()
    }
    # Seeds are a deterministic function of (sample size index, replicate), so
    # a run reproduces exactly regardless of how it is parallelized.
    base <- seed + (i - 1L) * n_sim * 7919L
    out <- lapply_maybe_parallel(seq_len(n_sim), function(j) {
      one_rep(n[i], j, base + j * 104729L)
    }, n_cores)
    reps <- c(reps, out)
  }

  rep_df <- do.call(rbind, lapply(reps, function(z) as.data.frame(z)))
  rownames(rep_df) <- NULL

  agg <- lapply(split(rep_df, rep_df$n), function(d) {
    p1 <- mean(d$detected)
    p2 <- mean(d$recovered_all)
    data.frame(
      n = d$n[1], n_sim = nrow(d),
      power = p1, se_power = sqrt(p1 * (1 - p1) / nrow(d)),
      power_recovery = p2, se_recovery = sqrt(p2 * (1 - p2) / nrow(d)),
      mean_dice = mean(d$mean_dice),
      mean_n_sig = mean(d$n_sig),
      false_subnet_rate = mean(d$false_subnet_rate, na.rm = TRUE))
  })
  power_df <- do.call(rbind, agg)
  power_df <- power_df[order(power_df$n), ]
  rownames(power_df) <- NULL

  by_cluster <- do.call(rbind, lapply(seq_len(n_clust), function(k) {
    col <- paste0("dice", k)
    do.call(rbind, lapply(split(rep_df, rep_df$n), function(d) {
      data.frame(n = d$n[1], cluster = k, size = cluster_size[k],
                 f2 = rep_len(f2, n_clust)[k],
                 power = mean(d[[col]] >= dice_cut),
                 mean_dice = mean(d[[col]]))
    }))
  }))
  by_cluster <- by_cluster[order(by_cluster$n, by_cluster$cluster), ]
  rownames(by_cluster) <- NULL

  structure(
    list(power = power_df, by_cluster = by_cluster, replicates = rep_df,
         params = list(n = n, n_nodes = n_nodes, cluster_size = cluster_size,
                       f2 = f2, rho_in = rho_in, rho_out = rho_out,
                       n_cov = n_cov, n_sim = n_sim, n_perm = n_perm,
                       alpha = alpha, dice_cut = dice_cut,
                       threshold_prob = threshold_prob, lambda = lambda,
                       objective = objective, tune_each = tune_each,
                       seed = seed)),
    class = "subnet_power")
}

#' Sample size required to reach a target power
#'
#' Interpolates a [power_curve()] to report the smallest sample size reaching a
#' target power.
#'
#' @param object A `subnet_power` object.
#' @param target Target power, e.g. `0.8`.
#' @param which Which power curve to read: `"recovery"` (default, the
#'   probability of recovering the planted subnetworks) or `"any"` (the
#'   probability of any detection).
#'
#' @return A single number: the interpolated sample size, `NA` with a warning
#'   if the curve never reaches `target` over the range simulated.
#'
#' @details
#' Linear interpolation between the two bracketing sample sizes. Extend the
#' grid rather than trusting an extrapolation if the target lies outside it.
#'
#' @examples
#' \donttest{
#' pw <- power_curve(n = c(40, 80, 120), n_nodes = 60, cluster_size = 12,
#'                   f2 = 0.2, n_sim = 20, n_perm = 99, seed = 1)
#' required_n(pw, target = 0.8)
#' }
#'
#' @export
required_n <- function(object, target = 0.8, which = c("recovery", "any")) {
  stopifnot(inherits(object, "subnet_power"))
  which <- match.arg(which)
  d <- object$power
  y <- if (which == "recovery") d$power_recovery else d$power
  x <- d$n

  hit <- which(y >= target)
  if (!length(hit)) {
    warning("Power does not reach ", target, " for n up to ", max(x),
            "; extend the sample-size grid.", call. = FALSE)
    return(NA_real_)
  }
  i <- min(hit)
  if (i == 1L) return(x[1])
  x0 <- x[i - 1L]; x1 <- x[i]
  y0 <- y[i - 1L]; y1 <- y[i]
  if (y1 == y0) return(x1)
  x0 + (target - y0) * (x1 - x0) / (y1 - y0)
}
