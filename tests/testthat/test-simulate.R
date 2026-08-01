test_that("simulate_fc returns a coherent object", {
  sim <- simulate_fc(n = 80, n_nodes = 50, cluster_size = c(12, 8),
                     f2 = c(0.2, 0.3), seed = 1)

  expect_s3_class(sim, "subnet_sim")
  expect_equal(dim(sim$W), c(50L, 50L))
  expect_equal(sim$W, t(sim$W))
  expect_true(all(diag(sim$W) == 0))
  expect_true(all(is.finite(vech(sim$W))))
  expect_equal(lengths(sim$truth$nodes), c(12L, 8L))
  expect_length(intersect(sim$truth$nodes[[1]], sim$truth$nodes[[2]]), 0L)
  expect_equal(sim$params$df, 80L - 2L)
})

test_that("truth nodes survive the node shuffle", {
  # The planted blocks must remain blocks after relabelling: within-block edges
  # should be far more often supra-threshold than the rest.
  sim <- simulate_fc(n = 150, n_nodes = 60, cluster_size = 15, f2 = 0.3,
                     rho_in = 1, rho_out = 0, seed = 2)
  nodes <- sim$truth$nodes[[1]]
  inside <- vech(sim$W[nodes, nodes])
  M <- matrix(0, 60, 60); M[nodes, nodes] <- 1; diag(M) <- 0
  outside <- vech(sim$W)[vech(M) == 0]

  expect_gt(mean(inside), 5 * mean(outside))
})

test_that("edge_signal is consistent with the planted blocks", {
  sim <- simulate_fc(n = 100, n_nodes = 40, cluster_size = 10, f2 = 0.2,
                     rho_in = 1, rho_out = 0, seed = 3)
  nodes <- sim$truth$nodes[[1]]
  M <- matrix(0, 40, 40); M[nodes, nodes] <- 1; diag(M) <- 0
  expect_equal(sim$truth$edge_signal, vech(M) == 1)
})

test_that("rho_in and rho_out set the right number of signal edges", {
  cs <- 20L
  sim <- simulate_fc(n = 100, n_nodes = 60, cluster_size = cs, f2 = 0.2,
                     rho_in = 0.5, rho_out = 0.1, seed = 4)
  m_in <- cs * (cs - 1) / 2
  m_out <- 60 * 59 / 2 - m_in
  expect_equal(sum(sim$truth$edge_signal),
               floor(0.5 * m_in) + floor(0.1 * m_out))
})

test_that("the fast path and the data path agree in distribution", {
  # Both must produce the same null behaviour: -log10 p uniform on the log
  # scale means p itself is uniform.
  pf <- replicate(40, {
    s <- simulate_fc(n = 60, n_nodes = 30, cluster_size = 5, f2 = 0,
                     rho_out = 0, method = "fast")
    mean(10^(-vech(s$W)) < 0.05)
  })
  pd <- replicate(40, {
    s <- simulate_fc(n = 60, n_nodes = 30, cluster_size = 5, f2 = 0,
                     rho_out = 0, method = "data")
    mean(10^(-vech(s$W)) < 0.05)
  })
  expect_equal(mean(pf), 0.05, tolerance = 0.02)
  expect_equal(mean(pd), 0.05, tolerance = 0.02)
})

test_that("effect size shows up at the intended magnitude", {
  # For a signal edge the expected t-statistic is sqrt(f2 * n) to first order,
  # so a bigger f2 must give systematically larger -log10 p.
  small <- simulate_fc(n = 200, n_nodes = 40, cluster_size = 15, f2 = 0.05,
                       rho_in = 1, rho_out = 0, seed = 5)
  large <- simulate_fc(n = 200, n_nodes = 40, cluster_size = 15, f2 = 0.30,
                       rho_in = 1, rho_out = 0, seed = 5)
  expect_gt(mean(vech(large$W)[large$truth$edge_signal]),
            mean(vech(small$W)[small$truth$edge_signal]))
})

test_that("simulate_fc validates its arguments", {
  expect_error(simulate_fc(n = 50, n_nodes = 20, cluster_size = 25),
               "exceeds `n_nodes`")
  expect_error(simulate_fc(n = 50, n_nodes = 20, cluster_size = 1),
               "at least 2")
  expect_error(simulate_fc(n = 4, n_nodes = 20, cluster_size = 5, n_cov = 5),
               "residual degrees of freedom")
})

test_that("seeding makes simulation reproducible", {
  a <- simulate_fc(n = 60, n_nodes = 30, cluster_size = 8, f2 = 0.2, seed = 7)
  b <- simulate_fc(n = 60, n_nodes = 30, cluster_size = 8, f2 = 0.2, seed = 7)
  expect_equal(a$W, b$W)
  expect_equal(a$truth$nodes, b$truth$nodes)
})
