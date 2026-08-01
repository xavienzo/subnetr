# Timing of the extraction core, the permutation test and a power analysis,
# across connectome sizes.
#
# Run with:  Rscript benchmarks/benchmark.R
#
# Cost is driven by the number of edges that survive screening rather than by
# the number of nodes, so the screening quantile matters more here than n_nodes
# does. Both quantiles below are reported for that reason.

suppressMessages(library(subnetr))

make_W <- function(n_nodes, size = max(5, n_nodes %/% 5), seed = 1) {
  set.seed(seed)
  W <- matrix(0, n_nodes, n_nodes)
  W[lower.tri(W)] <- rexp(n_nodes * (n_nodes - 1) / 2)
  W <- W + t(W)
  W[1:size, 1:size] <- W[1:size, 1:size] + 5
  diag(W) <- 0
  W
}

timeit <- function(expr, min_time = 0.5, max_reps = 2000L) {
  # `expr` must be re-evaluated on every iteration; a plain promise would be
  # evaluated once and then cached, silently reporting near-zero timings.
  e <- substitute(expr)
  env <- parent.frame()
  reps <- 0L
  t0 <- Sys.time()
  repeat {
    eval(e, env)
    reps <- reps + 1L
    el <- as.numeric(Sys.time() - t0, units = "secs")
    if (el >= min_time || reps >= max_reps) break
  }
  el / reps
}

cat(sprintf("subnetr benchmarks -- R %s, %s, OpenMP: %s\n\n",
            getRversion(), Sys.info()[["sysname"]], has_openmp()))

cat("Single extraction\n\n")
rows <- list()
for (N in c(50, 100, 200, 400, 800)) {
  W <- make_W(N)
  for (q in c(0.95, 0.99)) {
    thr <- as.numeric(quantile(vech(W), q))
    rows[[length(rows) + 1L]] <- data.frame(
      nodes = N,
      screen = q,
      edges_kept = sum(vech(W) >= thr),
      ms = round(1000 * timeit(subnet_extract(W, thr, lambda = 0.6)), 3))
  }
}
print(do.call(rbind, rows), row.names = FALSE)

cat("\n\nPermutation test, 200 permutations at fixed parameters\n\n")
rows <- list()
for (N in c(100, 200, 400, 800)) {
  W <- make_W(N)
  thr <- as.numeric(quantile(vech(W), 0.99))
  el <- system.time(
    subnet(W, threshold = thr, lambda = 0.6, n_perm = 200, seed = 1)
  )[["elapsed"]]
  rows[[length(rows) + 1L]] <- data.frame(
    nodes = N, seconds = round(el, 3), ms_per_perm = round(1000 * el / 200, 2))
}
print(do.call(rbind, rows), row.names = FALSE)

cat("\n\nFull analysis including threshold and lambda tuning\n\n")
rows <- list()
for (N in c(100, 200, 400)) {
  W <- make_W(N)
  el <- system.time(subnet(W, n_perm = 500, seed = 1))[["elapsed"]]
  rows[[length(rows) + 1L]] <- data.frame(nodes = N, seconds = round(el, 2))
}
print(do.call(rbind, rows), row.names = FALSE)

cat("\n\nPower analysis, 200 replicates x 199 permutations\n\n")
el <- system.time(
  power_curve(n = 100, n_nodes = 100, cluster_size = 20, f2 = 0.08,
              n_sim = 200, n_perm = 199, seed = 1, progress = FALSE)
)[["elapsed"]]
cat(sprintf("  %.1f s single-core (%d extractions)\n", el, 200 * 200))
