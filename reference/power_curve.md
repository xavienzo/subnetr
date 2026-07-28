# Power analysis for subnetwork detection

Estimates, by simulation, the probability that a study of a given size
detects a predictor-associated subnetwork of a given size and effect
size.

## Usage

``` r
power_curve(
  n,
  n_nodes = 100L,
  cluster_size = 20L,
  f2 = 0.05,
  rho_in = 0.9,
  rho_out = 0.02,
  f2_out = NULL,
  n_cov = 0L,
  n_sim = 200L,
  n_perm = 199L,
  alpha = 0.05,
  threshold_prob = NULL,
  lambda = 0.6,
  objective = c("generalized", "density"),
  min_size = 3L,
  max_clusters = 25L,
  tune_each = FALSE,
  tune = list(),
  dice_cut = 0.5,
  n_cores = 1L,
  seed = NULL,
  progress = TRUE
)
```

## Arguments

- n:

  Integer vector of sample sizes to evaluate.

- n_nodes, cluster_size, f2, rho_in, rho_out, f2_out, n_cov:

  Passed to
  [`simulate_fc()`](https://xavienzo.github.io/subnetr/reference/simulate_fc.md),
  describing the connectome and the planted effect.

- n_sim:

  Replicates per sample size.

- n_perm:

  Permutations per replicate. `199` gives p-value resolution `0.005`,
  ample for `alpha = 0.05`.

- alpha:

  Family-wise significance level.

- threshold_prob:

  Quantile of the observed edge weights used as the screening threshold
  when `tune_each = FALSE`. `NULL` (the default) scales it with the
  design so that screening retains roughly twice as many edges as the
  planted subnetworks contain, capped at the top 10%. See the section on
  the screening threshold – this is the single most consequential
  setting in a power analysis.

- lambda:

  Objective parameter used when `tune_each = FALSE`.

- objective, min_size, max_clusters:

  Passed to
  [`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md).

- tune_each:

  Retune the threshold and `lambda` in every replicate with
  [`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md),
  and pay the corresponding selection cost in the permutation null. This
  simulates the procedure a real analysis runs, which is the reason to
  use it; it is considerably slower, and the power it reports can land
  either side of a fixed threshold.

- tune:

  A named list of extra arguments for
  [`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md),
  used only when `tune_each = TRUE`. Shrinking its grid or its `n_perm`
  is the main way to make `tune_each` affordable.

- dice_cut:

  Minimum Sorensen-Dice overlap for a planted subnetwork to count as
  recovered.

- n_cores:

  Cores used to run replicates in parallel.

- seed:

  Optional integer seed. Results are reproducible and do not depend on
  `n_cores`.

- progress:

  Print a line per sample size.

## Value

An object of class `subnet_power`:

- power:

  Data frame with one row per sample size: `n`, `power`,
  `power_recovery`, their Monte Carlo standard errors `se_power` and
  `se_recovery`, `mean_dice`, `mean_n_sig` and `false_subnet_rate`, the
  average proportion of significant subnetworks that match no planted
  one.

- by_cluster:

  Data frame with one row per sample size and planted subnetwork, giving
  that subnetwork's own recovery probability and mean overlap.

- replicates:

  Per-replicate results, for diagnostics.

- params:

  The settings used.

## Details

For each sample size and each replicate the function simulates a
connectome study with
[`simulate_fc()`](https://xavienzo.github.io/subnetr/reference/simulate_fc.md),
runs the full
[`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md)
pipeline including its permutation test, and records what was found. Two
quantities are reported because they answer different questions:

- `power`:

  the probability of declaring *any* subnetwork significant. This is the
  conventional notion of power, but on its own it can be satisfied by a
  subnetwork that has little to do with the planted one.

- `power_recovery`:

  the probability that every planted subnetwork is matched by a
  significant one with Sorensen-Dice overlap at least `dice_cut`. This
  is the quantity you want when the scientific claim is about *which*
  nodes are involved, and it is always the more demanding of the two.

The gap between them is informative: when `power` is high but
`power_recovery` is not, the study is large enough to see that something
is there but not to localize it.

Replicates are independent, so the Monte Carlo standard error of each
proportion is reported and shrinks as `1 / sqrt(n_sim)`. Roughly 200
replicates give a standard error near 0.035 at 50% power.

## The screening threshold

Nothing else in the analysis moves power as much as `threshold_prob`. A
subnetwork of `c` nodes spans `c(c-1)/2` edges; if screening retains
fewer edges than that across the whole graph, the subnetwork cannot
survive as a block no matter how large the sample. In an 80-node
connectome with a 20-node planted subnetwork at f2 = 0.06 and n = 120,
recovery probability runs

|                  |          |
|------------------|----------|
| `threshold_prob` | recovery |
| 0.99             | 0.30     |
| 0.975            | 0.95     |
| 0.95             | 1.00     |

purely from where the threshold sits: the subnetwork spans 190 edges,
and at the 99th percentile only 32 edges survive screening in the entire
graph. The `NULL` default scales the threshold with the design to keep
this from silently dominating the answer, but a reported power is a
statement about a *particular* screening choice and should be quoted
alongside it.

If you intend to tune the threshold at analysis time rather than fix it
in advance, set `tune_each = TRUE` so the simulation runs the procedure
you will actually use. That is slower, and it can come out either higher
or lower than a fixed threshold – permutations now pay the same search
cost, which costs power, but the threshold adapts to each data set,
which gains it. The reason to use it is fidelity, not conservatism.

## See also

[`required_n()`](https://xavienzo.github.io/subnetr/reference/required_n.md),
[`plot.subnet_power()`](https://xavienzo.github.io/subnetr/reference/plot.subnet_power.md),
[`simulate_fc()`](https://xavienzo.github.io/subnetr/reference/simulate_fc.md),
[`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md)

## Examples

``` r
# \donttest{
pw <- power_curve(n = c(50, 100), n_nodes = 60, cluster_size = 12,
                  f2 = 0.2, n_sim = 20, n_perm = 99, seed = 1)
#> n = 50  (20 replicates)
#> n = 100  (20 replicates)
pw
#> Power analysis for subnetwork detection
#> 
#>   nodes           : 60
#>   planted subnets : 12 (sizes)
#>   effect size f2  : 0.2
#>   rho_in / rho_out: 0.90 / 0.020
#>   replicates      : 20   permutations: 99   alpha: 0.05
#>   Dice cutoff     : 0.50
#> 
#>    n         power      recovery  dice n_sig
#>   50 1.000 (0.000) 1.000 (0.000) 0.982     1
#>  100 1.000 (0.000) 1.000 (0.000) 1.000     1
#> 
#> Power and recovery show the estimate with its Monte Carlo standard
#> error in parentheses.
# }
```
