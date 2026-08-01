test_that("vech and unvech round-trip and match MATLAB squareform order", {
  W <- matrix(0, 5, 5)
  W[lower.tri(W)] <- 1:10
  W <- W + t(W)

  expect_equal(vech(W), as.numeric(1:10))
  expect_equal(unvech(vech(W)), W)

  # Column-major lower-triangle order: the first entry is (2,1), the second
  # (3,1), and so on -- the same order MATLAB's squareform produces.
  expect_equal(vech(W)[1], W[2, 1])
  expect_equal(vech(W)[2], W[3, 1])
  expect_equal(vech(W)[5], W[3, 2])
})

test_that("unvech rejects impossible lengths", {
  expect_error(unvech(1:5), "not a valid number of edges")
})

test_that("check_adjacency enforces its contract", {
  expect_error(subnet_extract(1:10, 1), "must be a numeric matrix")
  expect_error(subnet_extract(matrix(1, 2, 3), 1), "must be square")
  expect_error(subnet_extract(matrix(1, 2, 2), 1), "at least 3 nodes")

  A <- matrix(runif(25), 5, 5)
  expect_error(subnet_extract(A, 0.5), "must be symmetric")

  B <- A + t(A); diag(B) <- 0
  B[1, 1] <- NA
  expect_error(subnet_extract(B, 0.5), "missing values")
})

test_that("perm_jobs covers every permutation of every grid point exactly once", {
  # This runs both branches on every platform. The OpenMP branch is otherwise
  # unreachable on a machine without OpenMP, and a gap in it silently leaves
  # part of the null unfilled rather than raising an error.
  for (use_omp in c(TRUE, FALSE)) {
    for (ng in c(1L, 2L, 5L, 25L)) {
      for (n_perm in c(1L, 49L, 999L)) {
        for (n_cores in c(1L, 2L, 4L, 8L)) {
          jobs <- subnetr:::perm_jobs(ng, n_perm, n_cores, use_omp)
          info <- sprintf("omp=%s ng=%d n_perm=%d n_cores=%d",
                          use_omp, ng, n_perm, n_cores)
          for (g in seq_len(ng)) {
            mine <- Filter(function(j) j$g == g, jobs)
            covered <- unlist(lapply(mine, function(j) {
              j$ch$offset + seq_len(j$ch$n)
            }))
            expect_equal(sort(covered), seq_len(n_perm), info = info)
          }
          expect_setequal(vapply(jobs, `[[`, numeric(1), "g"), seq_len(ng))
        }
      }
    }
  }
})

test_that("perm_jobs issues one job per grid point when OpenMP threads inside", {
  # With OpenMP the C++ loop is already parallel, so splitting permutations
  # across R workers as well would oversubscribe the machine.
  jobs <- subnetr:::perm_jobs(ng = 1L, n_perm = 500L, n_cores = 8L,
                              use_omp = TRUE)
  expect_length(jobs, 1L)
  expect_equal(jobs[[1]]$ch$n, 500L)

  # Without it, the same work is spread over the available workers.
  jobs <- subnetr:::perm_jobs(ng = 1L, n_perm = 500L, n_cores = 8L,
                              use_omp = FALSE)
  expect_gt(length(jobs), 1L)
  expect_equal(sum(vapply(jobs, function(j) j$ch$n, numeric(1))), 500)
})

test_that("perm_chunks partitions permutations exactly once", {
  for (np in c(1L, 7L, 100L)) {
    for (nc in c(1L, 3L, 8L, 200L)) {
      ch <- subnetr:::perm_chunks(np, nc)
      expect_equal(sum(vapply(ch, `[[`, numeric(1), "n")), np)
      offs <- vapply(ch, `[[`, numeric(1), "offset")
      ns <- vapply(ch, `[[`, numeric(1), "n")
      covered <- unlist(Map(function(o, n) o + seq_len(n), offs, ns))
      expect_equal(sort(covered), seq_len(np))
    }
  }
})
