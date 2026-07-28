# Tune the screening threshold and objective parameter

Selects the pair (`threshold`, `lambda`) that yields the most structured
partition, over a grid of candidate values.

## Usage

``` r
tune_subnet(
  W,
  lambda_grid = seq(0.5, 0.9, by = 0.1),
  probs = c(0.9, 0.95, 0.975, 0.99, 0.995),
  threshold_grid = NULL,
  criterion = c("calibrated", "likelihood"),
  n_perm = 25L,
  objective = c("generalized", "density"),
  min_size = 3L,
  max_clusters = 25L,
  n_cores = 1L,
  seed = NULL
)
```

## Arguments

- W:

  Symmetric numeric matrix of edge weights.

- lambda_grid:

  Candidate values for the size-penalty exponent. Ignored when
  `objective` is not `"generalized"`, in which case only the first value
  is used.

- probs:

  Quantiles of the observed edge weights used to build the candidate
  thresholds.

- threshold_grid:

  Candidate thresholds. Overrides `probs` when given.

- criterion:

  `"calibrated"` (default) or `"likelihood"`; see Details.

- n_perm:

  Permutations per grid point for the calibrated criterion. 25 is
  usually enough to rank grid points reliably.

- objective, min_size, max_clusters:

  Passed to
  [`subnet_extract()`](https://xavienzo.github.io/subnetr/reference/subnet_extract.md).

- n_cores:

  Cores used to evaluate the grid.

- seed:

  Optional seed for the calibration permutations.

## Value

An object of class `subnet_tuning`:

- threshold, lambda:

  The selected values.

- grid:

  Data frame with one row per grid point holding `lambda`, `prob`,
  `threshold`, `n_clusters`, `lr` and (when calibrated) `lr_null`,
  `lr_null_sd` and `z`.

- criterion:

  The criterion used.

## Details

Every candidate pair is scored by the log-likelihood ratio of a block
model against a homogeneous model. Writing \\k_i\\ for the number of
supra-threshold edges inside extracted subnetwork \\i\\, \\m_i\\ for the
number of edges it could contain, and \\\pi\\ for the overall
supra-threshold rate, the fitted model gives each subnetwork its own
edge probability \\\hat\pi_i = k_i / m_i\\ and everything else a shared
\\\hat\pi_0\\. The statistic is the resulting gain in Bernoulli
log-likelihood over the single-\\\pi\\ model.

That statistic cannot be compared across thresholds as it stands,
because changing the threshold changes the binary data being modelled –
a higher threshold mechanically produces a sparser graph and a different
likelihood scale. `criterion = "calibrated"` (the default) puts every
grid point on a common scale by standardizing its likelihood ratio
against a null obtained by permuting the edge weights and rerunning the
same extraction: \$\$z(\lambda, r) = \frac{LR\_{obs}(\lambda, r) -
\mathrm{mean}(LR\_{null})}{\mathrm{sd}(LR\_{null})}.\$\$ The grid point
with the largest `z` is selected. All grid points share the same
permutations, so they are compared under common random numbers. This
costs `n_perm * length(lambda_grid) * length(probs)` extractions, which
the compiled peeling code makes cheap.

Set `criterion = "likelihood"` to fall back on the uncalibrated
statistic, which is faster but biased toward the extreme ends of the
threshold grid.

## See also

[`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md),
[`subnet_extract()`](https://xavienzo.github.io/subnetr/reference/subnet_extract.md)

## Examples

``` r
sim <- simulate_fc(n = 100, n_nodes = 60, cluster_size = 12, f2 = 0.2,
                   seed = 7)
tn <- tune_subnet(sim$W, n_perm = 10, seed = 1)
tn
#> Threshold / lambda tuning
#>   criterion : calibrated (10 permutations per grid point)
#>   grid      : 5 lambda x 5 threshold
#>   selected  : lambda = 0.70, threshold = 2.283
#> 
#> Top grid points:
#>  lambda prob threshold n_clusters    lr     z
#>     0.7 0.95     2.283          8 180.4 21.63
#>     0.6 0.95     2.283          8 180.4 17.37
#>     0.8 0.90     1.291         11 113.8 15.86
#>     0.7 0.90     1.291         10 153.1 15.47
#>     0.6 0.90     1.291          8 150.1 11.96
```
