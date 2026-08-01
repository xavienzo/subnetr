test_that("tune_subnet returns a usable grid and selection", {
  sim <- simulate_fc(n = 100, n_nodes = 50, cluster_size = 12, f2 = 0.25,
                     seed = 31)
  tn <- tune_subnet(sim$W, n_perm = 8, seed = 1)

  expect_s3_class(tn, "subnet_tuning")
  expect_equal(nrow(tn$grid), 5L * 5L)
  expect_true(tn$lambda %in% tn$grid$lambda)
  expect_true(tn$threshold %in% tn$grid$threshold)
  expect_true(all(tn$grid$lr >= 0))
  expect_true(any(is.finite(tn$grid$z)))
})

test_that("the likelihood criterion needs no permutations", {
  sim <- simulate_fc(n = 100, n_nodes = 50, cluster_size = 12, f2 = 0.25,
                     seed = 32)
  tn <- tune_subnet(sim$W, criterion = "likelihood")
  expect_identical(tn$criterion, "likelihood")
  expect_true(all(is.na(tn$grid$z)))
  expect_equal(tn$n_perm, 0L)
})

test_that("tuning finds a threshold that separates signal from noise", {
  sim <- simulate_fc(n = 150, n_nodes = 60, cluster_size = 15, f2 = 0.3,
                     seed = 33)
  tn <- tune_subnet(sim$W, n_perm = 10, seed = 2)
  part <- subnet_extract(sim$W, tn$threshold, tn$lambda)
  expect_gte(dice(part$blocks[[1]], sim$truth$nodes[[1]]), 0.8)
})

test_that("lambda_grid collapses for objectives that ignore lambda", {
  sim <- simulate_fc(n = 80, n_nodes = 40, cluster_size = 10, f2 = 0.2,
                     seed = 34)
  tn <- tune_subnet(sim$W, objective = "density", criterion = "likelihood")
  expect_equal(length(unique(tn$grid$lambda)), 1L)
})

test_that("tune_subnet validates its arguments", {
  sim <- simulate_fc(n = 60, n_nodes = 30, cluster_size = 8, f2 = 0.2, seed = 35)
  expect_error(tune_subnet(sim$W, probs = c(0.5, 1.5)), "strictly between")
  expect_error(tune_subnet(sim$W, lambda_grid = c(0.5, 2)), "lie in \\[0, 1\\]")
  expect_error(tune_subnet(sim$W, n_perm = 2), "at least 3")
})

test_that("partition_lr is zero for a homogeneous partition", {
  # One block covering every edge with the overall rate must give no gain.
  k <- matrix(50, 1, 1); m <- matrix(100, 1, 1)
  expect_equal(subnetr:::partition_lr(k, m, 0L, 50, 100), 0)
})

test_that("partition_lr grows as a block departs from the overall rate", {
  M <- 1000; K <- 100
  weak <- subnetr:::partition_lr(matrix(c(20, 0), 1), matrix(c(100, 0), 1),
                                 1L, K, M)
  strong <- subnetr:::partition_lr(matrix(c(80, 0), 1), matrix(c(100, 0), 1),
                                   1L, K, M)
  expect_gt(strong, weak)
  expect_gt(weak, 0)
})

test_that("power_curve produces sane output that rises with sample size", {
  skip_on_cran()
  # A weak effect and a wide sample-size range, so the increase in power is
  # far larger than the Monte Carlo noise at this number of replicates.
  pw <- power_curve(n = c(30, 200), n_nodes = 40, cluster_size = 10,
                    f2 = 0.04, n_sim = 40, n_perm = 99, seed = 41,
                    progress = FALSE)

  expect_s3_class(pw, "subnet_power")
  expect_equal(nrow(pw$power), 2L)
  expect_equal(nrow(pw$by_cluster), 2L)
  expect_true(all(pw$power$power >= 0 & pw$power$power <= 1))
  # Recovery is a strictly stronger requirement than any detection.
  expect_true(all(pw$power$power_recovery <= pw$power$power))
  expect_gt(pw$power$power_recovery[2], pw$power$power_recovery[1] + 0.05)
  expect_gt(pw$power$mean_dice[2], 2 * pw$power$mean_dice[1])
  expect_equal(nrow(pw$replicates), 80L)
})

test_that("power_curve handles several planted subnetworks", {
  skip_on_cran()
  pw <- power_curve(n = 100, n_nodes = 50, cluster_size = c(12, 8),
                    f2 = c(0.3, 0.3), n_sim = 20, n_perm = 49, seed = 42,
                    progress = FALSE)
  expect_equal(nrow(pw$by_cluster), 2L)
  expect_equal(pw$by_cluster$size, c(12L, 8L))
})

test_that("required_n interpolates and warns when the target is unreachable", {
  skip_on_cran()
  pw <- power_curve(n = c(30, 60, 120), n_nodes = 40, cluster_size = 10,
                    f2 = 0.15, n_sim = 30, n_perm = 99, seed = 43,
                    progress = FALSE)
  n80 <- required_n(pw, target = 0.8, which = "any")
  expect_true(is.na(n80) || (n80 >= 30 && n80 <= 120))
  expect_warning(required_n(pw, target = 1.01), "does not reach")
})
