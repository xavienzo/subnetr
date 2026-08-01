#' Edge-wise association statistics for a connectivity study
#'
#' Regresses every edge of a connectome on a predictor of interest while
#' adjusting for nuisance covariates, and returns the resulting evidence as a
#' weighted adjacency matrix suitable for [subnet()].
#'
#' @details
#' The model fitted at each edge \eqn{e} is
#' \deqn{y_e = \beta_{0e} + \beta_e x + Z \gamma_e + \varepsilon_e.}
#' Rather than looping over edges, the predictor and the connectivity matrix
#' are both residualized on `cbind(1, covariates)` once, after which the
#' Frisch-Waugh-Lovell theorem gives every edge's slope, standard error and
#' t-statistic in a handful of vectorized operations. The cost is one QR
#' factorization of an `n * (p + 1)` matrix plus two BLAS-level products,
#' independent of the number of edges.
#'
#' P-values are computed on the log scale, so evidence far beyond double
#' precision -- routine in connectome studies with thousands of edges -- is
#' represented exactly instead of saturating at `Inf`.
#'
#' @param fc Connectivity data. Either an `n_subject` by `n_edge` matrix whose
#'   rows are vectorized lower triangles (see [vech()]), or an
#'   `n_node` by `n_node` by `n_subject` array of connectivity matrices.
#' @param x Numeric predictor of interest, length `n_subject`.
#' @param covariates Optional numeric matrix or data frame of nuisance
#'   covariates with `n_subject` rows. An intercept is always included and must
#'   not be supplied.
#' @param value What to return in the matrix: `"nlog10p"` (default),
#'   `"t"`, or `"p"`.
#'
#' @return A symmetric `n_node` by `n_node` matrix with a zero diagonal,
#'   carrying attributes `df` (residual degrees of freedom), `t` (the vector of
#'   edge t-statistics) and `beta` (the vector of edge slopes).
#'
#' @seealso [subnet()], [simulate_fc()]
#'
#' @examples
#' set.seed(1)
#' n <- 60; p <- 20; m <- p * (p - 1) / 2
#' x <- rnorm(n)
#' fc <- matrix(rnorm(n * m), n, m)
#' fc[, 1:10] <- fc[, 1:10] + 0.6 * x        # plant a signal
#' W <- edge_stats(fc, x, covariates = cbind(age = rnorm(n)))
#' dim(W)
#' round(max(W), 2)
#'
#' @export
edge_stats <- function(fc, x, covariates = NULL,
                       value = c("nlog10p", "t", "p")) {
  value <- match.arg(value)

  if (is.array(fc) && length(dim(fc)) == 3L) {
    d <- dim(fc)
    if (d[1] != d[2]) stop("`fc` array must be n_node x n_node x n_subject.",
                           call. = FALSE)
    idx <- which(lower.tri(matrix(0, d[1], d[1])))
    Y <- t(matrix(fc, d[1] * d[2], d[3])[idx, , drop = FALSE])
  } else if (is.matrix(fc)) {
    Y <- fc
  } else {
    stop("`fc` must be a matrix or a 3-dimensional array.", call. = FALSE)
  }
  storage.mode(Y) <- "double"

  n <- nrow(Y)
  M <- ncol(Y)
  if (length(x) != n) {
    stop("`x` must have one value per subject (", n, ").", call. = FALSE)
  }
  if (anyNA(Y) || anyNA(x)) {
    stop("Missing values are not supported; drop or impute them first.",
         call. = FALSE)
  }

  Z0 <- if (is.null(covariates)) {
    matrix(1, n, 1)
  } else {
    cbind(1, data.matrix(covariates))
  }
  if (nrow(Z0) != n) {
    stop("`covariates` must have one row per subject (", n, ").", call. = FALSE)
  }

  qrZ <- qr(Z0)
  p <- qrZ$rank
  df <- n - p - 1L
  if (df < 1L) {
    stop("Not enough residual degrees of freedom (n = ", n,
         ", covariates = ", p, ").", call. = FALSE)
  }

  xr <- qr.resid(qrZ, as.numeric(x))
  Sxx <- sum(xr^2)
  if (Sxx <= .Machine$double.eps) {
    stop("`x` is collinear with the covariates.", call. = FALSE)
  }
  Yr <- qr.resid(qrZ, Y)

  beta <- as.numeric(crossprod(Yr, xr)) / Sxx
  rss <- colSums(Yr^2) - beta^2 * Sxx
  rss[rss < 0] <- 0
  se <- sqrt(rss / df / Sxx)

  tstat <- ifelse(se > 0, beta / se, 0)

  W <- unvech(switch(
    value,
    t = tstat,
    p = 2 * pt(-abs(tstat), df),
    nlog10p = -(log(2) + pt(-abs(tstat), df, log.p = TRUE)) / log(10)))

  nnode <- nrow(W)
  if (nnode * (nnode - 1) / 2 != M) {
    stop("`fc` has ", M, " columns, which is not a valid number of edges.",
         call. = FALSE)
  }

  attr(W, "df") <- df
  attr(W, "t") <- tstat
  attr(W, "beta") <- beta
  # Recorded so that plots can label their colour scale with the quantity
  # actually being shown rather than a generic "weight".
  attr(W, "statistic") <- switch(value,
    t = "t statistic",
    p = "p-value",
    nlog10p = expression(-log[10](p)))
  W
}
