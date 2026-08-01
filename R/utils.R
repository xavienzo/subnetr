## Internal helpers -----------------------------------------------------------

OBJ_CODES <- c(generalized = 0L, density = 1L)

match_objective <- function(objective) {
  objective <- match.arg(objective, names(OBJ_CODES))
  OBJ_CODES[[objective]]
}

#' Convert between a symmetric matrix and its vectorized lower triangle
#'
#' `vech()` extracts the strict lower triangle of a symmetric matrix in
#' column-major order and `unvech()` rebuilds the matrix. This is the ordering
#' the package uses for every edge vector, and it matches `W[lower.tri(W)]` in
#' base R and `squareform()` in MATLAB, so edge vectors can be exchanged with
#' either without reindexing.
#'
#' @param W A symmetric numeric matrix with a zero diagonal.
#' @param v A numeric vector of length `n * (n - 1) / 2`.
#'
#' @return `vech()` returns a numeric vector of length `n * (n - 1) / 2`;
#'   `unvech()` returns a symmetric numeric matrix with a zero diagonal.
#'
#' @examples
#' W <- matrix(0, 4, 4)
#' W[lower.tri(W)] <- 1:6
#' W <- W + t(W)
#' vech(W)
#' identical(unvech(vech(W)), W)
#'
#' @export
vech <- function(W) {
  W <- as.matrix(W)
  W[lower.tri(W)]
}

#' @rdname vech
#' @export
unvech <- function(v) {
  m <- length(v)
  n <- (1 + sqrt(1 + 8 * m)) / 2
  if (abs(n - round(n)) > 1e-8) {
    stop("`v` has length ", m, ", which is not a valid number of edges.",
         call. = FALSE)
  }
  n <- as.integer(round(n))
  W <- matrix(0, n, n)
  W[lower.tri(W)] <- v
  W + t(W)
}

check_adjacency <- function(W, arg = "W") {
  if (!is.matrix(W) || !is.numeric(W)) {
    stop("`", arg, "` must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(W) != ncol(W)) {
    stop("`", arg, "` must be square.", call. = FALSE)
  }
  if (nrow(W) < 3L) {
    stop("`", arg, "` must have at least 3 nodes.", call. = FALSE)
  }
  if (anyNA(W)) {
    stop("`", arg, "` must not contain missing values.", call. = FALSE)
  }
  asym <- max(abs(W - t(W)))
  if (asym > 1e-8 * max(1, max(abs(W)))) {
    stop("`", arg, "` must be symmetric (max asymmetry ",
         format(asym, digits = 3), ").", call. = FALSE)
  }
  lab <- attr(W, "statistic")
  storage.mode(W) <- "double"
  diag(W) <- 0
  attr(W, "statistic") <- lab
  W
}

## Upper-tail binomial log p-value for each block.
##
## A block of `size` nodes spans m = size*(size-1)/2 possible edges of which
## `k` exceed the screening threshold. Under exchangeability each edge is
## supra-threshold with probability `pi0`, so the evidence that the block is
## denser than chance is P(Binom(m, pi0) >= k), returned on the log scale
## because the values routinely underflow double precision.
block_logp <- function(k, m, pi0) {
  out <- numeric(length(k))
  pi0 <- min(max(pi0, .Machine$double.eps), 1 - .Machine$double.eps)
  ok <- m >= 1 & k >= 1
  if (any(ok)) {
    out[ok] <- pbinom(k[ok] - 1, size = m[ok], prob = pi0,
                      lower.tail = FALSE, log.p = TRUE)
  }
  out
}

## Split the (grid point, permutation) work of the permutation null into jobs.
##
## Every grid point must receive every permutation exactly once. The only
## question is how the permutations for a grid point are executed: as a single
## call whose C++ loop is threaded by OpenMP, or split across several forked R
## workers when OpenMP is unavailable. Getting this wrong drops permutations
## silently and corrupts the null, so the split is isolated here and tested
## directly rather than inferred from a run on one platform.
perm_jobs <- function(ng, n_perm, n_cores, use_omp) {
  jobs <- list()
  for (g in seq_len(ng)) {
    chunks <- if (use_omp || ng >= n_cores) {
      list(list(offset = 0L, n = n_perm))
    } else {
      perm_chunks(n_perm, ceiling(n_cores / ng))
    }
    for (ch in chunks) jobs[[length(jobs) + 1L]] <- list(g = g, ch = ch)
  }
  jobs
}

## Split `n_perm` permutations into `n_cores` contiguous chunks. Because each
## permutation seeds its own RNG stream from its absolute index, chunking never
## changes the result.
perm_chunks <- function(n_perm, n_cores) {
  n_cores <- max(1L, min(as.integer(n_cores), n_perm))
  bounds <- floor(seq(0, n_perm, length.out = n_cores + 1L))
  keep <- diff(bounds) > 0
  Map(function(from, size) list(offset = from, n = size),
      bounds[-length(bounds)][keep], diff(bounds)[keep])
}

## Run `f` over `x`, in parallel when asked for and available.
lapply_maybe_parallel <- function(x, f, n_cores) {
  if (n_cores <= 1L || !requireNamespace("parallel", quietly = TRUE) ||
      .Platform$OS.type == "windows") {
    return(lapply(x, f))
  }
  parallel::mclapply(x, f, mc.cores = n_cores)
}

#' Report whether the compiled code has OpenMP support
#'
#' The permutation loop is parallelized with OpenMP where the toolchain
#' provides it. When it does not -- most notably the default Apple clang on
#' macOS -- `subnet_test()` falls back to forked R workers, so the `n_cores`
#' argument still delivers a speedup.
#'
#' @return A single logical value.
#'
#' @examples
#' has_openmp()
#'
#' @export
has_openmp <- function() subnetr_has_openmp()
