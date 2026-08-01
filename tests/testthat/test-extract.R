make_planted <- function(n_nodes = 40, size = 10, seed = 1) {
  set.seed(seed)
  W <- matrix(0, n_nodes, n_nodes)
  W[lower.tri(W)] <- runif(n_nodes * (n_nodes - 1) / 2, 0, 1)
  W <- W + t(W)
  W[1:size, 1:size] <- W[1:size, 1:size] + 4
  diag(W) <- 0
  W
}

test_that("lambda trades subnetwork size against density", {
  W <- make_planted(n_nodes = 50, size = 12, seed = 9)
  thr <- quantile(vech(W), 0.85)

  # The size penalty is n^(2*lambda), so raising lambda must not enlarge the
  # extracted block, and the extremes must differ.
  sizes <- vapply(c(0.05, 0.3, 0.6, 0.95), function(l) {
    subnet_extract(W, thr, lambda = l)$sizes[1]
  }, numeric(1))
  expect_false(is.unsorted(rev(sizes)))
  expect_gt(sizes[1], sizes[length(sizes)])
})

test_that("lambda = 0.5 maximizes average degree", {
  # w / n^(2 * 0.5) = w / n, so the selected block must be the one maximizing
  # average weighted degree over the peeling sequence.
  set.seed(21)
  N <- 40L
  W <- matrix(0, N, N)
  W[lower.tri(W)] <- rexp(N * (N - 1) / 2)
  W <- W + t(W)
  W[1:10, 1:10] <- W[1:10, 1:10] + 4
  diag(W) <- 0
  thr <- quantile(vech(W), 0.9)

  got <- subnet_extract(W, thr, lambda = 0.5, min_size = 2L)$blocks[[1]]

  A <- W; A[A < thr] <- 0
  best <- NULL; best_score <- -Inf
  alive <- seq_len(N)
  repeat {
    sub <- A[alive, alive, drop = FALSE]
    n <- length(alive)
    if (n < 2L) break
    score <- sum(sub) / 2 / n
    if (score > best_score) { best_score <- score; best <- alive }
    alive <- alive[-which.min(colSums(sub))]
  }
  expect_setequal(got, best)
})

test_that("extraction recovers a strongly planted subnetwork", {
  size <- 12
  W <- make_planted(n_nodes = 50, size = size, seed = 3)
  part <- subnet_extract(W, threshold = quantile(vech(W), 0.9))
  expect_setequal(part$blocks[[1]], seq_len(size))
})

test_that("blocks partition the nodes exactly once", {
  W <- make_planted(n_nodes = 45, size = 11, seed = 4)
  part <- subnet_extract(W, threshold = quantile(vech(W), 0.9))
  expect_setequal(unlist(part$blocks), seq_len(45))
  expect_equal(sum(part$sizes), 45L)
  expect_equal(length(part$order), 45L)
  expect_false(anyDuplicated(part$order) > 0)
})

test_that("min_size is respected", {
  W <- make_planted(n_nodes = 45, size = 11, seed = 6)
  part <- subnet_extract(W, threshold = quantile(vech(W), 0.9), min_size = 8L)
  nk <- part$n_clusters
  if (nk > 0) expect_true(all(part$sizes[seq_len(nk)] >= 8L))
})

test_that("reported edge counts and densities are correct", {
  W <- make_planted(n_nodes = 40, size = 10, seed = 8)
  thr <- quantile(vech(W), 0.9)
  part <- subnet_extract(W, threshold = thr)

  for (i in seq_along(part$blocks)) {
    b <- part$blocks[[i]]
    if (length(b) < 2) next
    sub <- W[b, b, drop = FALSE]
    expect_equal(part$n_edge[i], sum(vech(sub) >= thr))
    expect_equal(part$density[i],
                 sum(vech(sub) >= thr) / (length(b) * (length(b) - 1) / 2))
  }
})

test_that("a graph with no supra-threshold edges yields no subnetworks", {
  W <- matrix(0, 20, 20)
  part <- subnet_extract(W, threshold = 1)
  expect_equal(part$n_clusters, 0L)
  expect_equal(part$sizes, 20L)
})

test_that("max_clusters caps the number of rounds", {
  W <- make_planted(n_nodes = 60, size = 10, seed = 12)
  part <- subnet_extract(W, threshold = quantile(vech(W), 0.8),
                         max_clusters = 2L)
  expect_lte(part$n_clusters, 2L)
  expect_setequal(unlist(part$blocks), seq_len(60))
})

test_that("extraction is invariant to node relabelling", {
  W <- make_planted(n_nodes = 45, size = 11, seed = 15)
  thr <- quantile(vech(W), 0.9)
  a <- subnet_extract(W, thr)

  set.seed(2)
  p <- sample.int(45)
  b <- subnet_extract(W[p, p], thr)
  expect_setequal(p[b$blocks[[1]]], a$blocks[[1]])
})
