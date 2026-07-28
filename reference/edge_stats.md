# Edge-wise association statistics for a connectivity study

Regresses every edge of a connectome on a predictor of interest while
adjusting for nuisance covariates, and returns the resulting evidence as
a weighted adjacency matrix suitable for
[`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md).

## Usage

``` r
edge_stats(fc, x, covariates = NULL, value = c("nlog10p", "t", "p"))
```

## Arguments

- fc:

  Connectivity data. Either an `n_subject` by `n_edge` matrix whose rows
  are vectorized lower triangles (see
  [`vech()`](https://xavienzo.github.io/subnetr/reference/vech.md)), or
  an `n_node` by `n_node` by `n_subject` array of connectivity matrices.

- x:

  Numeric predictor of interest, length `n_subject`.

- covariates:

  Optional numeric matrix or data frame of nuisance covariates with
  `n_subject` rows. An intercept is always included and must not be
  supplied.

- value:

  What to return in the matrix: `"nlog10p"` (default), `"t"`, or `"p"`.

## Value

A symmetric `n_node` by `n_node` matrix with a zero diagonal, carrying
attributes `df` (residual degrees of freedom), `t` (the vector of edge
t-statistics) and `beta` (the vector of edge slopes).

## Details

The model fitted at each edge \\e\\ is \$\$y_e = \beta\_{0e} + \beta_e
x + Z \gamma_e + \varepsilon_e.\$\$ Rather than looping over edges, the
predictor and the connectivity matrix are both residualized on
`cbind(1, covariates)` once, after which the Frisch-Waugh-Lovell theorem
gives every edge's slope, standard error and t-statistic in a handful of
vectorized operations. The cost is one QR factorization of an
`n * (p + 1)` matrix plus two BLAS-level products, independent of the
number of edges.

P-values are computed on the log scale, so evidence far beyond double
precision – routine in connectome studies with thousands of edges – is
represented exactly instead of saturating at `Inf`.

## See also

[`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md),
[`simulate_fc()`](https://xavienzo.github.io/subnetr/reference/simulate_fc.md)

## Examples

``` r
set.seed(1)
n <- 60; p <- 20; m <- p * (p - 1) / 2
x <- rnorm(n)
fc <- matrix(rnorm(n * m), n, m)
fc[, 1:10] <- fc[, 1:10] + 0.6 * x        # plant a signal
W <- edge_stats(fc, x, covariates = cbind(age = rnorm(n)))
dim(W)
#> [1] 20 20
round(max(W), 2)
#> [1] 6.36
```
