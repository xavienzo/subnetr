# Convert between a symmetric matrix and its vectorized lower triangle

`vech()` extracts the strict lower triangle of a symmetric matrix in
column-major order and `unvech()` rebuilds the matrix. This is the
ordering the package uses for every edge vector, and it matches
`W[lower.tri(W)]` in base R and `squareform()` in MATLAB, so edge
vectors can be exchanged with either without reindexing.

## Usage

``` r
vech(W)

unvech(v)
```

## Arguments

- W:

  A symmetric numeric matrix with a zero diagonal.

- v:

  A numeric vector of length `n * (n - 1) / 2`.

## Value

`vech()` returns a numeric vector of length `n * (n - 1) / 2`;
`unvech()` returns a symmetric numeric matrix with a zero diagonal.

## Examples

``` r
W <- matrix(0, 4, 4)
W[lower.tri(W)] <- 1:6
W <- W + t(W)
vech(W)
#> [1] 1 2 3 4 5 6
identical(unvech(vech(W)), W)
#> [1] TRUE
```
