test_that("subnet detects a strong planted subnetwork and localizes it", {
  sim <- simulate_fc(n = 120, n_nodes = 60, cluster_size = 12, f2 = 0.25,
                     seed = 21)
  fit <- subnet(sim$W, n_perm = 199, seed = 1)

  expect_s3_class(fit, "subnet")
  expect_true(any(fit$significant))
  expect_gte(dice(fit$subnetworks[[1]], sim$truth$nodes[[1]]), 0.8)
  expect_true(fit$significant[1])
})

test_that("p-values are bounded and monotone in the statistic", {
  sim <- simulate_fc(n = 100, n_nodes = 50, cluster_size = 12, f2 = 0.2,
                     seed = 22)
  fit <- subnet(sim$W, n_perm = 99, seed = 2)

  expect_true(all(fit$p_value >= 1 / (fit$n_perm + 1)))
  expect_true(all(fit$p_value <= 1))
  # A more extreme (more negative) statistic can never get a larger p-value.
  o <- order(fit$statistic)
  expect_false(is.unsorted(fit$p_value[o]))
})

test_that("results are reproducible and independent of n_cores", {
  sim <- simulate_fc(n = 100, n_nodes = 50, cluster_size = 12, f2 = 0.2,
                     seed = 23)
  a <- subnet(sim$W, threshold = quantile(vech(sim$W), 0.95), lambda = 0.6,
              n_perm = 199, seed = 99)
  b <- subnet(sim$W, threshold = quantile(vech(sim$W), 0.95), lambda = 0.6,
              n_perm = 199, seed = 99)
  expect_equal(a$p_value, b$p_value)
  expect_equal(a$null_statistic, b$null_statistic)

  skip_on_os("windows")
  skip_if_not_installed("parallel")
  cc <- subnet(sim$W, threshold = quantile(vech(sim$W), 0.95), lambda = 0.6,
               n_perm = 199, seed = 99, n_cores = 2)
  # Each permutation seeds its own stream from its absolute index, so chunking
  # across workers must not change a single value.
  expect_equal(cc$null_statistic, a$null_statistic)
  expect_equal(cc$p_value, a$p_value)
})

test_that("the null is calibrated when no signal is present", {
  skip_on_cran()
  nrep <- 200
  p <- vapply(seq_len(nrep), function(i) {
    sim <- simulate_fc(n = 80, n_nodes = 50, cluster_size = 10, f2 = 0,
                       rho_out = 0, seed = 500 + i)
    fit <- subnet(sim$W, threshold = quantile(vech(sim$W), 0.95),
                  lambda = 0.6, n_perm = 199, seed = 900 + i)
    min(fit$p_value)
  }, numeric(1))

  # Monte Carlo se at alpha = 0.05 with 200 replicates is about 0.015.
  expect_lt(abs(mean(p < 0.05) - 0.05), 0.05)
  expect_lt(abs(mean(p < 0.20) - 0.20), 0.09)
})

test_that("tuning on the same data inflates the uncorrected null", {
  skip_on_cran()
  # Documents the flaw that null = "retune" exists to fix: with the naive
  # null, selecting the threshold from the data roughly doubles type-I error.
  nrep <- 150
  run <- function(mode) {
    vapply(seq_len(nrep), function(i) {
      sim <- simulate_fc(n = 80, n_nodes = 50, cluster_size = 10, f2 = 0,
                         rho_out = 0, seed = 700 + i)
      fit <- subnet(sim$W, n_perm = 99, seed = 800 + i,
                    tune = list(n_perm = 10), null = mode)
      min(fit$p_value)
    }, numeric(1))
  }
  naive <- mean(run("selected") < 0.05)
  fixed <- mean(run("retune") < 0.05)
  expect_gt(naive, fixed)
  expect_lt(fixed, 0.12)
})

test_that("null modes are wired to the right permutation grid", {
  sim <- simulate_fc(n = 100, n_nodes = 40, cluster_size = 10, f2 = 0.2,
                     seed = 24)
  fit <- subnet(sim$W, n_perm = 49, seed = 3, tune = list(n_perm = 5))
  expect_identical(fit$null, "retune")
  expect_gt(nrow(fit$null_grid), 1L)

  fixed <- subnet(sim$W, threshold = 3, lambda = 0.6, n_perm = 49, seed = 3)
  expect_identical(fixed$null, "selected")
  expect_equal(nrow(fixed$null_grid), 1L)

  expect_error(subnet(sim$W, threshold = 3, lambda = 0.6, null = "retune",
                      n_perm = 49), "requires tuning")
})

test_that("reusing a tuning object matches tuning inside subnet()", {
  sim <- simulate_fc(n = 100, n_nodes = 40, cluster_size = 10, f2 = 0.2,
                     seed = 26)
  tn <- tune_subnet(sim$W, n_perm = 5, seed = 11)

  # Reusing the object must give exactly the analysis that tuning internally
  # would have given, including the grid the null searches. The internal call
  # needs the same tuning seed, or it draws its own and may land on a
  # different grid point.
  reused <- subnet(sim$W, tuning = tn, n_perm = 49, seed = 7)
  inline <- subnet(sim$W, n_perm = 49, seed = 7,
                   tune = list(n_perm = 5, seed = 11))

  expect_identical(reused$null, "retune")
  expect_equal(reused$threshold, tn$threshold)
  expect_equal(reused$lambda, tn$lambda)
  expect_equal(nrow(reused$null_grid), nrow(tn$grid))
  expect_equal(reused$p_value, inline$p_value)
})

test_that("passing tuned values as numbers drops the correction", {
  # Documents why `tuning =` exists: the same numbers supplied plainly carry
  # no record of having been selected, so the null is built as if they had
  # been fixed in advance.
  sim <- simulate_fc(n = 100, n_nodes = 40, cluster_size = 10, f2 = 0.2,
                     seed = 27)
  tn <- tune_subnet(sim$W, n_perm = 5, seed = 12)

  plain <- subnet(sim$W, threshold = tn$threshold, lambda = tn$lambda,
                  n_perm = 49, seed = 7)
  expect_identical(plain$null, "selected")
  expect_equal(nrow(plain$null_grid), 1L)
})

test_that("subnet() validates the tuning argument", {
  sim <- simulate_fc(n = 80, n_nodes = 40, cluster_size = 10, f2 = 0.2,
                     seed = 28)
  tn <- tune_subnet(sim$W, n_perm = 5, seed = 13)

  expect_error(subnet(sim$W, tuning = tn, threshold = 2, n_perm = 9),
               "not both")
  expect_error(subnet(sim$W, tuning = list(a = 1), n_perm = 9),
               "must be a `subnet_tuning` object")
  expect_error(subnet(sim$W, tuning = tn, min_size = 8L, n_perm = 9),
               "different extraction settings")
  expect_error(subnet(sim$W, tuning = tn, max_clusters = 5L, n_perm = 9),
               "different extraction settings")
})

test_that("membership and as.data.frame agree with the fit", {
  sim <- simulate_fc(n = 120, n_nodes = 50, cluster_size = 12, f2 = 0.25,
                     seed = 25)
  fit <- subnet(sim$W, n_perm = 99, seed = 4)

  df <- as.data.frame(fit)
  expect_equal(nrow(df), length(fit$size))
  expect_equal(df$size, fit$size)

  mb <- membership(fit)
  expect_length(mb, 50L)
  for (i in which(fit$significant)) {
    expect_setequal(which(mb == i), fit$subnetworks[[i]])
  }
  expect_equal(sum(mb > 0), sum(fit$size[fit$significant]))
})

test_that("dice behaves at the boundaries", {
  expect_equal(dice(1:10, 1:10), 1)
  expect_equal(dice(1:10, 11:20), 0)
  expect_equal(dice(integer(0), 1:5), 0)
  expect_equal(dice(1:10, 6:15), 0.5)
})
