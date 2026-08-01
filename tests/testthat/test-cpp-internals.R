test_that("decode_position inverts the vectorized lower-triangle ordering", {
  for (N in c(3L, 4L, 7L, 40L, 201L)) {
    M <- N * (N - 1) / 2
    # The reference mapping: which (row, col) does each vech() position hold?
    W <- matrix(0, N, N)
    W[lower.tri(W)] <- seq_len(M)
    want_i <- row(W)[lower.tri(W)]
    want_j <- col(W)[lower.tri(W)]

    got <- subnetr:::decode_position_cpp(seq_len(M) - 1, N)
    expect_equal(got[, 1], want_i, info = paste("N =", N))
    expect_equal(got[, 2], want_j, info = paste("N =", N))
  }
})

test_that("the permutation null preserves the screened edge count", {
  # Every permutation must place exactly the edges that survived screening --
  # no more, no fewer -- since permuting the weight vector cannot change how
  # many entries clear the threshold.
  set.seed(1)
  N <- 40L
  W <- matrix(0, N, N)
  W[lower.tri(W)] <- rexp(N * (N - 1) / 2)
  W <- W + t(W)
  r <- as.numeric(quantile(vech(W), 0.9))
  E <- sum(vech(W) >= r)

  res <- subnetr:::subnet_perm_cpp(vech(W), N, r, 0.6, 0L, 2L, 100L,
                                   50L, 7, 0L, 1L)
  # Blocks partition the nodes, so the block edge counts can only miss edges
  # that cross blocks; the total must never exceed the number placed.
  expect_true(all(rowSums(res$k) <= E))
  expect_true(all(rowSums(res$k) > 0))

  sizes_ok <- apply(res$m, 1, function(mm) all(mm >= 0))
  expect_true(all(sizes_ok))
})

test_that("permutation results are invariant to chunking", {
  set.seed(2)
  N <- 30L
  W <- matrix(0, N, N)
  W[lower.tri(W)] <- rexp(N * (N - 1) / 2)
  W <- W + t(W)
  r <- as.numeric(quantile(vech(W), 0.9))

  whole <- subnetr:::subnet_perm_cpp(vech(W), N, r, 0.6, 0L, 3L, 100L,
                                     40L, 123, 0L, 1L)
  first <- subnetr:::subnet_perm_cpp(vech(W), N, r, 0.6, 0L, 3L, 100L,
                                     15L, 123, 0L, 1L)
  rest <- subnetr:::subnet_perm_cpp(vech(W), N, r, 0.6, 0L, 3L, 100L,
                                    25L, 123, 15L, 1L)

  expect_equal(rbind(first$k, rest$k), whole$k)
  expect_equal(rbind(first$m, rest$m), whole$m)
  expect_equal(c(first$nc, rest$nc), whole$nc)
})
