# Detect predictor-associated subnetworks with a permutation test

The main entry point. Screens a weighted connectivity matrix at a
threshold, extracts candidate subnetworks by greedy peeling, and assigns
each of them a family-wise-error-controlled p-value from a permutation
null.

## Usage

``` r
subnet(
  W,
  threshold = NULL,
  lambda = NULL,
  alpha = 0.05,
  n_perm = 1000L,
  objective = c("generalized", "density"),
  min_size = 3L,
  max_clusters = 25L,
  tune = list(),
  null = c("auto", "retune", "grid", "selected"),
  n_cores = 1L,
  seed = NULL,
  verbose = FALSE
)
```

## Arguments

- W:

  Symmetric numeric matrix of edge weights, typically \\-\log\_{10}\\
  p-values from
  [`edge_stats()`](https://xavienzo.github.io/subnetr/reference/edge_stats.md).

- threshold:

  Screening threshold. If `NULL` (default) it is chosen by
  [`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md).

- lambda:

  Size-penalty exponent. If `NULL` (default) it is chosen by
  [`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md).

- alpha:

  Family-wise significance level.

- n_perm:

  Number of permutations. The resolution of a p-value is
  `1 / (n_perm + 1)`, so `n_perm` should comfortably exceed `1 / alpha`.

- objective, min_size, max_clusters:

  Passed to
  [`subnet_extract()`](https://xavienzo.github.io/subnetr/reference/subnet_extract.md).

- tune:

  A named list of extra arguments for
  [`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md),
  used only when `threshold` or `lambda` is `NULL`.

- null:

  How to build the permutation null: `"auto"` (default) maximizes over
  the tuning grid when tuning was performed and over the selected
  parameters otherwise, `"grid"` always maximizes over the grid, and
  `"selected"` always uses the selected parameters only. See the section
  on tuning and validity.

- n_cores:

  Number of cores for the permutation loop. Uses OpenMP threads when the
  package was built with OpenMP support (see
  [`has_openmp()`](https://xavienzo.github.io/subnetr/reference/has_openmp.md))
  and forked R workers otherwise.

- seed:

  Optional integer seed for the permutation stream. When `NULL` the seed
  is drawn from R's RNG, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) makes a run
  reproducible. Results never depend on `n_cores`.

- verbose:

  Print progress messages.

## Value

An object of class `subnet`:

- subnetworks:

  List of integer node-index vectors, densest first.

- size, n_edge, density:

  Per-subnetwork node count, supra-threshold edge count, and
  within-block edge density.

- statistic:

  Per-subnetwork log binomial tail probability. More negative means
  denser than chance.

- p_value:

  Family-wise-error-controlled permutation p-value.

- significant:

  Logical, `p_value < alpha`.

- null_statistic:

  The `n_perm` permutation maxima.

- partition:

  The full
  [`subnet_extract()`](https://xavienzo.github.io/subnetr/reference/subnet_extract.md)
  partition, background block included.

- tuning:

  The
  [`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md)
  result, or `NULL` if both parameters were supplied.

- threshold, lambda, objective, alpha, n_perm, W:

  Settings and input.

## Details

Each block of the partition is scored by the upper-tail binomial
probability of seeing at least as many supra-threshold edges inside it
as were observed, if edges were placed at random with the graph's
overall supra-threshold rate. Because a block's size enters through the
binomial sample size, blocks of very different sizes are placed on a
common scale, which is what makes the next step legitimate.

The null is built by repeatedly permuting the edge weights across the
graph – destroying any subnetwork structure while preserving the
marginal distribution of edge weights – rerunning the *entire*
extraction, and recording the single most extreme block statistic
produced. Comparing each observed block against this distribution of
maxima is a Westfall-Young max-statistic procedure, so the reported
p-values control the family-wise error rate across all blocks the
algorithm returns, with no further multiplicity correction. It also
accounts for the selection effect of having chosen the blocks by
optimizing density in the first place: the null statistic is produced by
the same greedy search, applied to data with no signal in it.

P-values use the add-one estimator \\(1 + \\\\T\_{null} \le T\_{obs}\\)
/ (M + 1)\\, which is never zero and keeps the test exact at level
`alpha`.

## Tuning and validity

When `threshold` and `lambda` are chosen from the same matrix that is
then tested, a null built at the chosen values alone is optimistic: the
search had the opportunity to land on whichever setting made the data
look most structured, and the null never gets that opportunity. In
simulations of the global null this roughly doubles the type-I error,
taking a nominal 5% test to about 11%.

`subnet()` closes the gap by letting every permutation run the same
parameter search the observed data ran. Under the default
`null = "auto"`, when tuning was performed each permuted data set is
scored at every grid point, the tuning rule picks that permutation's own
`(lambda, threshold)`, and the maximum is taken over the blocks found
there. This mirrors the analysis step for step, so the test is
calibrated rather than merely valid. All grid points reuse the same
permuted data sets, so every selection happens within a permutation.

The reason this costs no more than one extraction pass per grid point is
that the tuning rule is a fixed function once its null moments are
known, and those moments describe the very distribution a permuted data
set is drawn from – so they are estimated once and reused, not
re-estimated inside each permutation.

Alternatives: `null = "grid"` takes the maximum over all grid points
instead, which is valid whatever the selection rule but noticeably
conservative (about 1.5% actual error for a nominal 5% test);
`null = "selected"` applies no correction at all and is appropriate only
when `threshold` and `lambda` were fixed in advance rather than tuned.

## See also

[`subnet_extract()`](https://xavienzo.github.io/subnetr/reference/subnet_extract.md),
[`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md),
[`edge_stats()`](https://xavienzo.github.io/subnetr/reference/edge_stats.md),
[`power_curve()`](https://xavienzo.github.io/subnetr/reference/power_curve.md),
[`plot.subnet()`](https://xavienzo.github.io/subnetr/reference/plot.subnet.md)

## Examples

``` r
sim <- simulate_fc(n = 120, n_nodes = 60, cluster_size = 12, f2 = 0.2,
                   seed = 42)
fit <- subnet(sim$W, n_perm = 199, seed = 1)
fit
#> Predictor-associated subnetworks
#> 
#>   nodes       : 60   edges: 1770
#>   objective   : generalized (lambda = 0.50)
#>   threshold   : 3.387   (5.0% of edges retained)
#>   permutations: 199   alpha: 0.05
#>   tuning      : calibrated criterion
#> 
#>  subnet size edges density p_value significant log10_stat
#>       1   12    55   0.833   0.005        TRUE      -59.6
#>       2    4     2   0.333   1.000       FALSE       -1.5
#>       3    8     5   0.179   1.000       FALSE       -1.9
#>       4   13     9   0.115   1.000       FALSE       -1.8
#>       5    3     1   0.333   1.000       FALSE       -0.8
#> 
#> 1 of 5 subnetworks significant at FWER 0.05.
#> P-values are max-statistic permutation p-values and are already
#> family-wise-error corrected; do not adjust them again.

# Recovery of the planted subnetwork
dice(fit$subnetworks[[1]], sim$truth$nodes[[1]])
#> [1] 1
```
