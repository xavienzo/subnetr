# Simulate a connectome study with planted predictor-associated subnetworks

Generates a weighted adjacency matrix of \\-\log\_{10}\\ p-values in
which one or more subnetworks carry an association with a predictor of
interest, together with the ground truth needed to score a detection
method.

## Usage

``` r
simulate_fc(
  n,
  n_nodes = 100L,
  cluster_size = 20L,
  f2 = 0.05,
  rho_in = 0.9,
  rho_out = 0.02,
  f2_out = NULL,
  n_cov = 0L,
  method = c("fast", "data"),
  edge_corr = 0,
  shuffle_nodes = TRUE,
  seed = NULL
)
```

## Arguments

- n:

  Number of subjects.

- n_nodes:

  Number of nodes in the connectome.

- cluster_size:

  Integer vector of planted subnetwork sizes.

- f2:

  Cohen's \\f^2\\ effect size for the affected edges, recycled across
  subnetworks so different subnetworks can carry different effects.

- rho_in:

  Proportion of within-subnetwork edges that carry the effect.

- rho_out:

  Proportion of outside edges that carry the effect.

- f2_out:

  Effect size for the scattered outside edges. Defaults to `f2[1]`.

- n_cov:

  Number of nuisance covariates adjusted for.

- method:

  `"fast"` or `"data"`; see the section on simulation method.

- edge_corr:

  Cross-edge noise correlation induced by a shared subject-level factor.
  `method = "data"` only.

- shuffle_nodes:

  Randomly relabel nodes so the planted subnetworks are scattered
  through the matrix.

- seed:

  Optional integer seed.

## Value

A list of class `subnet_sim`:

- W:

  `n_nodes` by `n_nodes` matrix of \\-\log\_{10}\\ p-values.

- truth:

  List with `nodes` (node indices of each planted subnetwork, in the
  returned labelling), `edge_signal` (logical over vectorized edges),
  and `labels` (block index per node, 0 for background).

- data:

  Subject-level `fc`, `x` and `covariates`, or `NULL` under
  `method = "fast"`.

- params:

  The settings used, including the residual `df`.

## Details

Nodes 1 through `sum(cluster_size)` are assigned to consecutive planted
subnetworks. A proportion `rho_in` of the edges inside each subnetwork
actually carry the effect (the rest are within-subnetwork nulls,
reflecting that a real subnetwork is not uniformly affected), and a
proportion `rho_out` of the edges outside every subnetwork carry it too
(scattered associations that are not part of any subnetwork). Node
labels are then shuffled so the planted structure is not sitting in the
corner of the matrix.

## Simulation method

`method = "fast"` (default) draws the edge t-statistics directly from
their exact sampling distribution rather than materializing
subject-level data. Under the linear model with fixed design, an edge
with Cohen's \\f^2\\ effect size has \$\$t \sim t\_{df}(\mathrm{ncp}),
\quad \mathrm{ncp} = \sqrt{f^2 S}, \quad S \sim \chi^2\_{n - 1 - q},\$\$
where \\df = n - 2 - q\\ and \\q\\ is the number of nuisance covariates.
`S` is drawn once per data set, which is what correctly couples all
edges through the shared predictor. This is exact, not an approximation,
and it avoids allocating an `n` by `n_edge` matrix – the reason power
analysis over thousands of replicates is practical.

`method = "data"` builds the subject-level edge matrix explicitly and
runs
[`edge_stats()`](https://xavienzo.github.io/subnetr/reference/edge_stats.md)
on it. Use it when you need the raw data, or when `edge_corr > 0` to
study how the test behaves under cross-edge dependence, which the fast
path deliberately excludes.

## See also

[`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md),
[`power_curve()`](https://xavienzo.github.io/subnetr/reference/power_curve.md),
[`edge_stats()`](https://xavienzo.github.io/subnetr/reference/edge_stats.md)

## Examples

``` r
sim <- simulate_fc(n = 100, n_nodes = 80, cluster_size = c(20, 10),
                   f2 = c(0.15, 0.25), seed = 3)
dim(sim$W)
#> [1] 80 80
lengths(sim$truth$nodes)
#> [1] 20 10
```
