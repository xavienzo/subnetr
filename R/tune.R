#' Tune the screening threshold and objective parameter
#'
#' Selects the pair (`threshold`, `lambda`) that yields the most structured
#' partition, over a grid of candidate values.
#'
#' @details
#' Every candidate pair is scored by the log-likelihood ratio of a block model
#' against a homogeneous model. Writing \eqn{k_i} for the number of
#' supra-threshold edges inside extracted subnetwork \eqn{i}, \eqn{m_i} for the
#' number of edges it could contain, and \eqn{\pi} for the overall
#' supra-threshold rate, the fitted model gives each subnetwork its own edge
#' probability \eqn{\hat\pi_i = k_i / m_i} and everything else a shared
#' \eqn{\hat\pi_0}. The statistic is the resulting gain in Bernoulli
#' log-likelihood over the single-\eqn{\pi} model.
#'
#' That statistic cannot be compared across thresholds as it stands, because
#' changing the threshold changes the binary data being modelled -- a higher
#' threshold mechanically produces a sparser graph and a different likelihood
#' scale. `criterion = "calibrated"` (the
#' default) puts every grid point on a common scale by standardizing its
#' likelihood ratio against a null obtained by permuting the edge weights and
#' rerunning the same extraction:
#' \deqn{z(\lambda, r) = \frac{LR_{obs}(\lambda, r) -
#'   \mathrm{mean}(LR_{null})}{\mathrm{sd}(LR_{null})}.}
#' The grid point with the largest `z` is selected. All grid points share the
#' same permutations, so they are compared under common random numbers. This
#' costs `n_perm * length(lambda_grid) * length(probs)` extractions, which the
#' compiled peeling code makes cheap.
#'
#' Set `criterion = "likelihood"` to fall back on the uncalibrated statistic,
#' which is faster but biased toward the extreme ends of the threshold grid.
#'
#' @param W Symmetric numeric matrix of edge weights.
#' @param lambda_grid Candidate values for the size-penalty exponent. Ignored
#'   when `objective` is not `"generalized"`, in which case only the first
#'   value is used.
#' @param probs Quantiles of the observed edge weights used to build the
#'   candidate thresholds.
#' @param threshold_grid Candidate thresholds. Overrides `probs` when given.
#' @param criterion `"calibrated"` (default) or `"likelihood"`; see Details.
#' @param n_perm Permutations per grid point for the calibrated criterion.
#'   25 is usually enough to rank grid points reliably.
#' @param objective,min_size,max_clusters Passed to [subnet_extract()].
#' @param n_cores Cores used to evaluate the grid.
#' @param seed Optional seed for the calibration permutations.
#'
#' @return An object of class `subnet_tuning`:
#'   \item{threshold, lambda}{The selected values.}
#'   \item{grid}{Data frame with one row per grid point holding `lambda`,
#'     `prob`, `threshold`, `n_clusters`, `lr` and (when calibrated) `lr_null`,
#'     `lr_null_sd` and `z`.}
#'   \item{criterion}{The criterion used.}
#'
#' @seealso [subnet()], [subnet_extract()]
#'
#' @examples
#' sim <- simulate_fc(n = 100, n_nodes = 60, cluster_size = 12, f2 = 0.2,
#'                    seed = 7)
#' tn <- tune_subnet(sim$W, n_perm = 10, seed = 1)
#' tn
#'
#' @export
tune_subnet <- function(W,
                        lambda_grid = seq(0.5, 0.9, by = 0.1),
                        probs = c(0.90, 0.95, 0.975, 0.99, 0.995),
                        threshold_grid = NULL,
                        criterion = c("calibrated", "likelihood"),
                        n_perm = 25L,
                        objective = c("generalized", "density"),
                        min_size = 3L, max_clusters = 25L,
                        n_cores = 1L, seed = NULL) {
  W <- check_adjacency(W)
  criterion <- match.arg(criterion)
  obj <- match_objective(objective)
  objective <- names(OBJ_CODES)[OBJ_CODES == obj]
  min_size <- max(2L, as.integer(min_size))
  max_clusters <- max(1L, as.integer(max_clusters))
  n_perm <- as.integer(n_perm)

  if (objective != "generalized") lambda_grid <- lambda_grid[1]
  if (any(lambda_grid < 0 | lambda_grid > 1)) {
    stop("`lambda_grid` values must lie in [0, 1].", call. = FALSE)
  }

  wv <- vech(W)
  N <- nrow(W)
  M_total <- length(wv)

  if (is.null(threshold_grid)) {
    if (any(probs <= 0 | probs >= 1)) {
      stop("`probs` values must lie strictly between 0 and 1.", call. = FALSE)
    }
    threshold_grid <- as.numeric(quantile(wv, probs, names = FALSE))
  } else {
    probs <- rep(NA_real_, length(threshold_grid))
  }
  keep <- !duplicated(threshold_grid)
  threshold_grid <- threshold_grid[keep]
  probs <- probs[keep]

  if (criterion == "calibrated" && n_perm < 3L) {
    stop("`n_perm` must be at least 3 for the calibrated criterion.",
         call. = FALSE)
  }
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  seed <- as.double(seed)

  grid <- expand.grid(lambda = lambda_grid, ti = seq_along(threshold_grid),
                      KEEP.OUT.ATTRS = FALSE)
  grid$threshold <- threshold_grid[grid$ti]
  grid$prob <- probs[grid$ti]

  res <- lapply_maybe_parallel(seq_len(nrow(grid)), function(g) {
    lam <- grid$lambda[g]
    r <- grid$threshold[g]
    K_total <- sum(wv >= r)

    fit <- subnet_extract_cpp(W, r, lam, obj, min_size, max_clusters)
    m_obs <- fit$sizes * (fit$sizes - 1) / 2
    lr_obs <- partition_lr(matrix(fit$n_edge, nrow = 1),
                           matrix(m_obs, nrow = 1),
                           fit$n_clusters, K_total, M_total)

    out <- list(n_clusters = fit$n_clusters, lr = lr_obs,
                lr_null = NA_real_, lr_null_sd = NA_real_, z = NA_real_)

    if (criterion == "calibrated") {
      # Same seed at every grid point: common random numbers.
      pn <- subnet_perm_cpp(wv, N, r, lam, obj, min_size, max_clusters,
                            n_perm, seed, 0L, 1L)
      lr_null <- partition_lr(pn$k, pn$m, pn$nc, K_total, M_total)
      mu <- mean(lr_null)
      s <- sd(lr_null)
      out$lr_null <- mu
      out$lr_null_sd <- s
      out$z <- if (is.finite(s) && s > 0) (lr_obs - mu) / s else NA_real_
    }
    out
  }, n_cores)

  grid$n_clusters <- vapply(res, `[[`, numeric(1), "n_clusters")
  grid$lr <- vapply(res, `[[`, numeric(1), "lr")
  grid$lr_null <- vapply(res, `[[`, numeric(1), "lr_null")
  grid$lr_null_sd <- vapply(res, `[[`, numeric(1), "lr_null_sd")
  grid$z <- vapply(res, `[[`, numeric(1), "z")
  grid$ti <- NULL

  if (criterion == "calibrated" && any(is.finite(grid$z))) {
    best <- which.max(replace(grid$z, !is.finite(grid$z), -Inf))
  } else {
    # Marginalize over thresholds to pick lambda, then pick the threshold.
    by_lambda <- tapply(grid$lr, grid$lambda, mean)
    lam_best <- as.numeric(names(by_lambda)[which.max(by_lambda)])
    cand <- which(grid$lambda == lam_best)
    best <- cand[which.max(grid$lr[cand])]
  }

  structure(
    list(threshold = grid$threshold[best],
         lambda = grid$lambda[best],
         grid = grid,
         best = best,
         criterion = criterion,
         objective = objective,
         # Recorded so that subnet() can verify a tuning result was produced
         # under the same extraction settings before reusing its grid.
         min_size = min_size,
         max_clusters = max_clusters,
         n_perm = if (criterion == "calibrated") n_perm else 0L),
    class = "subnet_tuning")
}

## Log-likelihood ratio of the block model against a homogeneous model, one
## value per row of `k`/`m`. `nc` gives the number of genuine subnetworks in
## each row; remaining columns are background or padding.
partition_lr <- function(k, m, nc, K_total, M_total) {
  k <- as.matrix(k); m <- as.matrix(m)
  nc <- rep_len(as.integer(nc), nrow(k))
  mask <- col(m) <= nc

  k_cl <- k * mask
  m_cl <- m * mask
  K1 <- rowSums(k_cl)
  M1 <- rowSums(m_cl)
  K0 <- K_total - K1
  M0 <- M_total - M1

  pi_ref <- K_total / M_total
  if (!is.finite(pi_ref) || pi_ref <= 0 || pi_ref >= 1) return(rep(0, nrow(k)))

  within <- rowSums(matrix(bern_lr(as.vector(k_cl), as.vector(m_cl), pi_ref),
                           nrow = nrow(k)))
  between <- bern_lr(K0, M0, pi_ref)
  out <- within + between
  out[!is.finite(out)] <- 0
  out
}

## k * log(pi_hat / pi_ref) + (m - k) * log((1 - pi_hat) / (1 - pi_ref)) with
## pi_hat = k / m, using the 0 * log(0) = 0 convention throughout.
bern_lr <- function(k, m, pi_ref) {
  out <- numeric(length(k))
  pos <- which(m > 0)
  if (!length(pos)) return(out)
  kk <- k[pos]; mm <- m[pos]
  a <- numeric(length(pos))
  i <- kk > 0
  a[i] <- kk[i] * (log(kk[i]) - log(mm[i]) - log(pi_ref))
  b <- numeric(length(pos))
  j <- (mm - kk) > 0
  b[j] <- (mm[j] - kk[j]) * (log(mm[j] - kk[j]) - log(mm[j]) - log1p(-pi_ref))
  out[pos] <- a + b
  out
}
