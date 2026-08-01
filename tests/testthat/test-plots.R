fit_for_plot <- function() {
  sim <- simulate_fc(n = 120, n_nodes = 40, cluster_size = 10, f2 = 0.25,
                     seed = 51)
  subnet(sim$W, threshold = quantile(vech(sim$W), 0.9), lambda = 0.6,
         n_perm = 49, seed = 1)
}

test_that("every matrix plot mode draws without error", {
  fit <- fit_for_plot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  expect_silent(plot(fit))
  expect_silent(plot(fit, what = "observed"))
  expect_silent(plot(fit, what = "both"))
  expect_silent(plot(fit, what = "both", legend = FALSE))
  expect_silent(plot(fit, significant_only = FALSE))
})

test_that("plotting restores the graphics state it found", {
  # The colour scale and the two-panel mode both need layout(); leaving either
  # in place would silently break whatever the user plots next.
  fit <- fit_for_plot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  for (mode in c("reordered", "observed", "both")) {
    before <- par(c("mfrow", "mar", "mfg"))
    plot(fit, what = mode)
    expect_equal(par(c("mfrow", "mar", "mfg")), before, info = mode)
    # A subsequent plot must occupy the whole device, not a stale sub-panel.
    expect_silent(plot(1:3, 1:3))
  }
})

test_that("the colour scale label follows the statistic", {
  fit <- fit_for_plot()
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  # simulate_fc() produces -log10 p and records it.
  expect_identical(attr(fit$W, "statistic"), expression(-log[10](p)))
  expect_silent(plot(fit, weight_label = "t statistic"))

  # A matrix carrying no label still plots, falling back to a generic one.
  bare <- fit
  attr(bare$W, "statistic") <- NULL
  expect_silent(plot(bare))
})

test_that("edge_stats labels each statistic it can return", {
  set.seed(52)
  n <- 40; m <- 10 * 9 / 2
  Y <- matrix(rnorm(n * m), n, m)
  x <- rnorm(n)
  expect_identical(attr(edge_stats(Y, x, value = "t"), "statistic"),
                   "t statistic")
  expect_identical(attr(edge_stats(Y, x, value = "p"), "statistic"),
                   "p-value")
  expect_identical(attr(edge_stats(Y, x, value = "nlog10p"), "statistic"),
                   expression(-log[10](p)))
})

test_that("power curves draw both series without warnings", {
  skip_on_cran()
  pw <- power_curve(n = c(40, 120), n_nodes = 40, cluster_size = 10,
                    f2 = 0.05, n_sim = 20, n_perm = 49, seed = 53,
                    progress = FALSE)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  # Power pinned at exactly 0 or 1 has zero standard error, which must not
  # produce a zero-length-arrow warning.
  expect_silent(plot(pw))
  expect_silent(plot(pw, which = "any"))
  expect_silent(plot(pw, which = "recovery", ci = FALSE))
  expect_silent(plot(pw, target = NULL))
})
