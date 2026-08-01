test_that("edge_stats matches lm() edge by edge", {
  set.seed(11)
  n <- 50; p <- 12; m <- p * (p - 1) / 2
  x <- rnorm(n)
  Z <- cbind(a = rnorm(n), b = rnorm(n))
  Y <- matrix(rnorm(n * m), n, m)
  Y[, 1:5] <- Y[, 1:5] + 0.8 * x

  W <- edge_stats(Y, x, covariates = Z, value = "t")
  got <- attr(W, "t")

  want <- vapply(seq_len(m), function(j) {
    summary(lm(Y[, j] ~ x + Z))$coefficients["x", "t value"]
  }, numeric(1))

  expect_equal(got, want, tolerance = 1e-10)
  expect_equal(attr(W, "df"), n - 3L - 1L)
})

test_that("edge_stats works without covariates", {
  set.seed(12)
  n <- 40; p <- 10; m <- p * (p - 1) / 2
  x <- rnorm(n)
  Y <- matrix(rnorm(n * m), n, m)

  W <- edge_stats(Y, x, value = "t")
  want <- vapply(seq_len(m), function(j) {
    summary(lm(Y[, j] ~ x))$coefficients["x", "t value"]
  }, numeric(1))
  expect_equal(attr(W, "t"), want, tolerance = 1e-10)
  expect_equal(attr(W, "df"), n - 2L)
})

test_that("the three value types agree with one another", {
  set.seed(13)
  n <- 40; p <- 10; m <- p * (p - 1) / 2
  x <- rnorm(n); Y <- matrix(rnorm(n * m), n, m)

  Wt <- edge_stats(Y, x, value = "t")
  Wp <- edge_stats(Y, x, value = "p")
  Wn <- edge_stats(Y, x, value = "nlog10p")
  df <- attr(Wt, "df")

  expect_equal(vech(Wp), 2 * pt(-abs(vech(Wt)), df))
  expect_equal(vech(Wn), -log10(vech(Wp)), tolerance = 1e-8)
})

test_that("log-scale p-values survive extreme evidence", {
  set.seed(14)
  n <- 200; p <- 8; m <- p * (p - 1) / 2
  x <- rnorm(n)
  Y <- matrix(rnorm(n * m, sd = 0.001), n, m) + 5 * x

  W <- edge_stats(Y, x)
  # -log10(p) far past the point where p itself underflows to zero.
  expect_true(all(is.finite(vech(W))))
  expect_gt(max(vech(W)), 100)
})

test_that("edge_stats accepts a 3-d array and agrees with the matrix form", {
  set.seed(15)
  n <- 30; p <- 8; m <- p * (p - 1) / 2
  x <- rnorm(n)
  Y <- matrix(rnorm(n * m), n, m)

  arr <- array(0, c(p, p, n))
  for (i in seq_len(n)) arr[, , i] <- unvech(Y[i, ])

  expect_equal(edge_stats(arr, x), edge_stats(Y, x))
})

test_that("edge_stats rejects bad input", {
  set.seed(16)
  n <- 20; m <- 10 * 9 / 2
  Y <- matrix(rnorm(n * m), n, m)
  expect_error(edge_stats(Y, rnorm(n - 1)), "one value per subject")
  expect_error(edge_stats(Y, rnorm(n), covariates = matrix(rnorm(n * 2), n)[-1, ]),
               "one row per subject")

  x <- rnorm(n)
  expect_error(edge_stats(Y, x, covariates = cbind(x)), "collinear")

  Yna <- Y; Yna[1, 1] <- NA
  expect_error(edge_stats(Yna, x), "Missing values")
})
