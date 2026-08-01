#' Detect predictor-associated subnetworks with a permutation test
#'
#' The main entry point. Screens a weighted connectivity matrix at a threshold,
#' extracts candidate subnetworks by greedy peeling, and assigns each of them a
#' family-wise-error-controlled p-value from a permutation null.
#'
#' @details
#' Each block of the partition is scored by the upper-tail binomial probability
#' of seeing at least as many supra-threshold edges inside it as were observed,
#' if edges were placed at random with the graph's overall supra-threshold rate.
#' Because a block's size enters through the binomial sample size, blocks of
#' very different sizes are placed on a common scale, which is what makes the
#' next step legitimate.
#'
#' The null is built by repeatedly permuting the edge weights across the graph
#' -- destroying any subnetwork structure while preserving the marginal
#' distribution of edge weights -- rerunning the *entire* extraction, and
#' recording the single most extreme block statistic produced. Comparing each
#' observed block against this distribution of maxima is a Westfall-Young
#' max-statistic procedure, so the reported p-values control the family-wise
#' error rate across all blocks the algorithm returns, with no further
#' multiplicity correction. It also accounts for the selection effect of having
#' chosen the blocks by optimizing density in the first place: the null
#' statistic is produced by the same greedy search, applied to data with no
#' signal in it.
#'
#' P-values use the add-one estimator \eqn{(1 + \#\{T_{null} \le T_{obs}\}) /
#' (M + 1)}, which is never zero and keeps the test exact at level `alpha`.
#'
#' @section Tuning and validity:
#' When `threshold` and `lambda` are chosen from the same matrix that is then
#' tested, a null built at the chosen values alone is optimistic: the search
#' had the opportunity to land on whichever setting made the data look most
#' structured, and the null never gets that opportunity. In simulations of the
#' global null this roughly doubles the type-I error, taking a nominal 5% test
#' to about 11%.
#'
#' `subnet()` closes the gap by letting every permutation run the same
#' parameter search the observed data ran. Under the default `null = "auto"`,
#' when tuning was performed each permuted data set is scored at every grid
#' point, the tuning rule picks that permutation's own `(lambda, threshold)`,
#' and the maximum is taken over the blocks found there. This mirrors the
#' analysis step for step, so the test is calibrated rather than merely valid.
#' All grid points reuse the same permuted data sets, so every selection
#' happens within a permutation.
#'
#' The reason this costs no more than one extraction pass per grid point is
#' that the tuning rule is a fixed function once its null moments are known,
#' and those moments describe the very distribution a permuted data set is
#' drawn from -- so they are estimated once and reused, not re-estimated inside
#' each permutation.
#'
#' Alternatives: `null = "grid"` takes the maximum over all grid points
#' instead, which is valid whatever the selection rule but noticeably
#' conservative (about 1.5% actual error for a nominal 5% test);
#' `null = "selected"` applies no correction at all and is appropriate only
#' when `threshold` and `lambda` were fixed in advance rather than tuned.
#'
#' Because the correction needs the grid that was searched, a tuning result
#' must be handed over whole. Running [tune_subnet()] and then passing its
#' selected values as plain numbers loses that record and silently reverts to
#' `null = "selected"`; pass `tuning =` instead.
#'
#' \preformatted{
#' tn <- tune_subnet(W)
#' subnet(W, threshold = tn$threshold, lambda = tn$lambda)  # uncorrected
#' subnet(W, tuning = tn)                                   # correct
#' }
#'
#' @param W Symmetric numeric matrix of edge weights, typically
#'   \eqn{-\log_{10}} p-values from [edge_stats()].
#' @param threshold Screening threshold. If `NULL` (default) it is chosen by
#'   [tune_subnet()].
#' @param lambda Size-penalty exponent. If `NULL` (default) it is chosen by
#'   [tune_subnet()].
#' @param alpha Family-wise significance level.
#' @param n_perm Number of permutations. The resolution of a p-value is
#'   `1 / (n_perm + 1)`, so `n_perm` should comfortably exceed `1 / alpha`.
#' @param objective,min_size,max_clusters Passed to [subnet_extract()].
#' @param tune A named list of extra arguments for [tune_subnet()], used only
#'   when `threshold` or `lambda` is `NULL`.
#' @param tuning A `subnet_tuning` object from a previous [tune_subnet()] call,
#'   to reuse rather than repeat the search. Cannot be combined with
#'   `threshold` or `lambda`. Prefer this to passing a tuning result's selected
#'   values as numbers: those carry no record of having been selected, so the
#'   null would be built as though the parameters had been fixed in advance.
#' @param null How to build the permutation null: `"auto"` (default) maximizes
#'   over the tuning grid when tuning was performed and over the selected
#'   parameters otherwise, `"grid"` always maximizes over the grid, and
#'   `"selected"` always uses the selected parameters only. See the section on
#'   tuning and validity.
#' @param n_cores Number of cores for the permutation loop. Uses OpenMP threads
#'   when the package was built with OpenMP support (see [has_openmp()]) and
#'   forked R workers otherwise.
#' @param seed Optional integer seed for the permutation stream. When `NULL`
#'   the seed is drawn from R's RNG, so `set.seed()` makes a run reproducible.
#'   Results never depend on `n_cores`.
#' @param verbose Print progress messages.
#'
#' @return An object of class `subnet`:
#'   \item{subnetworks}{List of integer node-index vectors, densest first.}
#'   \item{size, n_edge, density}{Per-subnetwork node count, supra-threshold
#'     edge count, and within-block edge density.}
#'   \item{statistic}{Per-subnetwork log binomial tail probability. More
#'     negative means denser than chance.}
#'   \item{p_value}{Family-wise-error-controlled permutation p-value.}
#'   \item{significant}{Logical, `p_value < alpha`.}
#'   \item{null_statistic}{The `n_perm` permutation maxima.}
#'   \item{partition}{The full [subnet_extract()] partition, background
#'     block included.}
#'   \item{tuning}{The [tune_subnet()] result, or `NULL` if both parameters
#'     were supplied.}
#'   \item{threshold, lambda, objective, alpha, n_perm, W}{Settings and input.}
#'
#' @seealso [subnet_extract()], [tune_subnet()], [edge_stats()],
#'   [power_curve()], [plot.subnet()]
#'
#' @examples
#' sim <- simulate_fc(n = 120, n_nodes = 60, cluster_size = 12, f2 = 0.2,
#'                    seed = 42)
#' fit <- subnet(sim$W, n_perm = 199, seed = 1)
#' fit
#'
#' # Recovery of the planted subnetwork
#' dice(fit$subnetworks[[1]], sim$truth$nodes[[1]])
#'
#' @export
subnet <- function(W, threshold = NULL, lambda = NULL, alpha = 0.05,
                   n_perm = 1000L,
                   objective = c("generalized", "density"),
                   min_size = 3L, max_clusters = 25L, tune = list(),
                   tuning = NULL,
                   null = c("auto", "retune", "grid", "selected"),
                   n_cores = 1L, seed = NULL, verbose = FALSE) {
  W <- check_adjacency(W)
  obj <- match_objective(objective)
  objective <- names(OBJ_CODES)[OBJ_CODES == obj]
  n_perm <- as.integer(n_perm)
  if (n_perm < 1L) stop("`n_perm` must be at least 1.", call. = FALSE)
  min_size <- max(2L, as.integer(min_size))
  max_clusters <- max(1L, as.integer(max_clusters))

  if (!is.null(tuning)) {
    # A tuning result carries the grid the search ranged over, which is what
    # the permutation null needs in order to repeat that search. Passing the
    # selected values as plain numbers instead would discard it.
    if (!inherits(tuning, "subnet_tuning")) {
      stop("`tuning` must be a `subnet_tuning` object from `tune_subnet()`.",
           call. = FALSE)
    }
    if (!is.null(threshold) || !is.null(lambda)) {
      stop("Supply either `tuning` or `threshold`/`lambda`, not both.",
           call. = FALSE)
    }
    if (!identical(tuning$objective, objective) ||
        !identical(as.integer(tuning$min_size), min_size) ||
        !identical(as.integer(tuning$max_clusters), max_clusters)) {
      stop("`tuning` was produced under different extraction settings ",
           "(objective, min_size or max_clusters); retune with the settings ",
           "used here.", call. = FALSE)
    }
    threshold <- tuning$threshold
    lambda <- tuning$lambda
  } else if (is.null(threshold) || is.null(lambda)) {
    if (verbose) message("Tuning threshold and lambda ...")
    tuning <- do.call(tune_subnet, c(
      list(W = W, objective = objective, min_size = min_size,
           max_clusters = max_clusters, n_cores = n_cores),
      tune))
    if (is.null(threshold)) threshold <- tuning$threshold
    if (is.null(lambda)) lambda <- tuning$lambda
  }

  if (verbose) {
    message(sprintf("Extracting subnetworks (threshold = %.3f, lambda = %.2f) ...",
                    threshold, lambda))
  }
  part <- subnet_extract(W, threshold = threshold, lambda = lambda,
                         objective = objective, min_size = min_size,
                         max_clusters = max_clusters)

  obs_logp <- partition_logp(part)

  # When the parameters were tuned on this same matrix, the null must be free
  # to search the same grid; otherwise the test inherits the optimism of the
  # selection. See the section on tuning and validity.
  null <- match.arg(null)
  if (null == "auto") null <- if (is.null(tuning)) "selected" else "retune"
  if (null != "selected" && is.null(tuning)) {
    stop("`null = \"", null, "\"` requires tuning; leave `threshold` or ",
         "`lambda` as NULL, or use `null = \"selected\"`.", call. = FALSE)
  }
  null_grid <- if (null == "selected") {
    data.frame(lambda = lambda, threshold = threshold)
  } else {
    tuning$grid[, c("lambda", "threshold")]
  }

  if (verbose) {
    message(sprintf("Running %d permutations over %d parameter setting%s ...",
                    n_perm, nrow(null_grid),
                    if (nrow(null_grid) == 1L) "" else "s"))
  }
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  null_stat <- perm_null(W = W, grid = null_grid, obj = obj,
                         min_size = min_size, max_clusters = max_clusters,
                         n_perm = n_perm, seed = seed, n_cores = n_cores,
                         mode = null, tuning = tuning)

  nk <- part$n_clusters
  idx <- seq_len(nk)
  pval <- vapply(idx, function(i) {
    (1 + sum(null_stat <= obs_logp[i])) / (n_perm + 1)
  }, numeric(1))

  structure(
    list(subnetworks = part$blocks[idx],
         size = part$sizes[idx],
         n_edge = part$n_edge[idx],
         density = part$density[idx],
         statistic = obs_logp[idx],
         p_value = pval,
         significant = pval < alpha,
         null_statistic = null_stat,
         partition = part,
         tuning = tuning,
         threshold = threshold,
         lambda = lambda,
         objective = objective,
         min_size = min_size,
         alpha = alpha,
         n_perm = n_perm,
         null = null,
         null_grid = null_grid,
         pi0 = part$pi0,
         W = W),
    class = "subnet")
}

## Permutation null: the distribution of the most extreme block statistic.
##
## `grid` is a data frame of (lambda, threshold) pairs and `mode` decides how
## the parameter search is accounted for:
##
##   "selected"  a single grid point, the ordinary max-over-blocks null. Correct
##               only when the parameters were not chosen from this matrix.
##   "grid"      maximum over blocks and over every grid point. Valid whatever
##               the selection rule was, but conservative, because the observed
##               statistic comes from one grid point while the null gets to
##               search all of them.
##   "retune"    each permutation runs the same selection rule the observed
##               data ran, picks its own grid point, and reports the maximum
##               over that point's blocks only. This mirrors the analysis
##               exactly and is therefore calibrated rather than conservative.
##
## The "retune" mode is affordable because the selection rule is a fixed
## function once tuning has estimated its null moments: those moments describe
## the null distribution of the likelihood ratio at each grid point, which is
## the same distribution a permuted data set is drawn from, so they can be
## reused rather than re-estimated inside every permutation.
##
## All grid points reuse the same permutations (identical seed), so maxima and
## selections happen within a permuted data set, never across unrelated ones.
perm_null <- function(W, grid, obj, min_size, max_clusters, n_perm, seed,
                      n_cores, mode = "selected", tuning = NULL) {
  wv <- vech(W)
  N <- nrow(W)
  M_total <- length(wv)
  n_cores <- max(1L, as.integer(n_cores))
  seed <- as.double(seed)
  ng <- nrow(grid)

  use_omp <- has_openmp()
  jobs <- perm_jobs(ng, n_perm, n_cores, use_omp)

  run <- function(job) {
    r <- grid$threshold[job$g]
    lam <- grid$lambda[job$g]
    res <- subnet_perm_cpp(wv, N, r, lam, obj, min_size, max_clusters,
                           job$ch$n, seed, job$ch$offset,
                           if (use_omp) n_cores else 1L)
    K_total <- sum(wv >= r)
    pi0 <- K_total / M_total
    lp <- matrix(block_logp(as.vector(res$k), as.vector(res$m), pi0),
                 nrow = nrow(res$k))
    list(g = job$g, offset = job$ch$offset,
         best = do.call(pmin, c(as.data.frame(lp), list(na.rm = TRUE))),
         lr = if (mode == "retune") {
           partition_lr(res$k, res$m, res$nc, K_total, M_total)
         } else NULL)
  }

  # With OpenMP the threading happens inside C++, so the jobs are run in
  # sequence here; otherwise they are spread over forked R workers.
  out <- if (use_omp) lapply(jobs, run) else
    lapply_maybe_parallel(jobs, run, n_cores)
  if (any(vapply(out, inherits, logical(1), "try-error"))) {
    stop("Permutation workers failed.", call. = FALSE)
  }

  best <- matrix(Inf, n_perm, ng)
  lr <- matrix(NA_real_, n_perm, ng)
  for (o in out) {
    idx <- o$offset + seq_along(o$best)
    best[idx, o$g] <- pmin(best[idx, o$g], o$best)
    if (!is.null(o$lr)) lr[idx, o$g] <- o$lr
  }

  if (mode == "retune") {
    pick <- select_grid_point(lr, grid, tuning$criterion,
                              tuning$grid$lr_null, tuning$grid$lr_null_sd)
    return(best[cbind(seq_len(n_perm), pick)])
  }
  if (ng == 1L) return(best[, 1L])
  do.call(pmin, c(as.data.frame(best), list(na.rm = TRUE)))
}

## Apply the tuning selection rule to each row of a matrix of per-grid-point
## likelihood ratios, returning the chosen column. Used identically on the
## observed data and inside permutations, which is what makes "retune" exact.
select_grid_point <- function(lr, grid, criterion, mu, sdv) {
  if (identical(criterion, "calibrated")) {
    z <- sweep(sweep(lr, 2, mu, "-"), 2, sdv, "/")
    z[!is.finite(z)] <- -Inf
    return(max.col(z, ties.method = "first"))
  }
  # "likelihood": average over thresholds to choose lambda, then choose the
  # threshold at that lambda.
  lam <- grid$lambda
  lev <- unique(lam)
  by_lam <- vapply(lev, function(l) rowMeans(lr[, lam == l, drop = FALSE]),
                   numeric(nrow(lr)))
  by_lam <- matrix(by_lam, nrow = nrow(lr))
  best_lam <- lev[max.col(by_lam, ties.method = "first")]
  vapply(seq_len(nrow(lr)), function(i) {
    cand <- which(lam == best_lam[i])
    cand[which.max(lr[i, cand])]
  }, integer(1))
}

#' Sorensen-Dice overlap between two node sets
#'
#' @param a,b Integer vectors of node indices.
#'
#' @return A number in \[0, 1\]; 1 means the two sets are identical and 0 that
#'   they are disjoint. Returns 0 when either set is empty.
#'
#' @examples
#' dice(1:20, 11:30)
#'
#' @export
dice <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0L || length(b) == 0L) return(0)
  2 * length(intersect(a, b)) / (length(a) + length(b))
}
