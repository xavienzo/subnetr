# Extract candidate subnetworks by greedy peeling

Partitions the nodes of a weighted graph into a sequence of increasingly
less dense blocks. Each round repeatedly deletes the minimum-degree
node, scores the set that remains after every deletion, and keeps the
best-scoring set as one subnetwork; the deleted nodes are then fed back
in for the next round. Everything left once the remainder carries no
supra-threshold weight is returned as a single background block.

## Usage

``` r
subnet_extract(
  W,
  threshold,
  lambda = 0.6,
  objective = c("sicers", "avg_degree", "density"),
  min_size = 3L,
  max_clusters = 25L
)
```

## Arguments

- W:

  Symmetric numeric matrix of edge weights, typically \\-\log\_{10}\\
  p-values from
  [`edge_stats()`](https://xavienzo.github.io/subnetr/reference/edge_stats.md).
  The diagonal is ignored.

- threshold:

  Screening threshold on the edge weights. Edges below it are set to
  zero before peeling.

- lambda:

  Objective tuning parameter in \[0, 1\]; see the section on choosing an
  objective. Ignored when `objective = "avg_degree"` or `"density"`.

- objective:

  One of `"sicers"` (default), `"avg_degree"` or `"density"`.

- min_size:

  Smallest number of nodes a block may contain. Blocks of two nodes are
  single edges and carry almost no evidence, so the default is 3.

- max_clusters:

  Cap on the number of extraction rounds. Peeling returns the densest
  block first and each subsequent round works on sparser leftovers, so
  the blocks that matter appear early; the cap simply stops the
  algorithm from grinding through a long tail of noise blocks. Raising
  it costs time and, in practice, does not change which blocks reach
  significance.

## Value

An object of class `subnet_partition`, a list with components

- blocks:

  List of integer vectors of node indices, densest block first. The
  final element is the background block and may be empty.

- sizes:

  Number of nodes in each block.

- n_edge:

  Number of supra-threshold edges inside each block.

- density:

  Proportion of possible within-block edges that are supra-threshold.

- n_clusters:

  Number of extracted subnetworks, excluding the background block.

- order:

  The node ordering that produces the block-diagonal arrangement, i.e.
  `blocks` concatenated.

- threshold, lambda, objective, min_size:

  The settings used.

- pi0:

  Proportion of all edges that are supra-threshold.

## Details

This is the extraction step only. It reports how dense each block is but
attaches no significance to it – use
[`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md) for
the full analysis with a calibrated permutation test.

## Choosing an objective

With `n` nodes spanning `m = n * (n - 1) / 2` possible edges and total
supra-threshold weight `w`, the objectives are

- `"sicers"`:

  \\(w/m)^\lambda \\ w^{1-\lambda}\\, the generalized density of Chen et
  al. (2023). `lambda = 1` targets pure density and therefore small
  tight cliques; `lambda = 0` targets total weight and therefore large
  diffuse blocks. Values in \[0.5, 0.9\] trade the two off.

- `"avg_degree"`:

  \\w/n\\, Charikar's classical densest-subgraph objective. This is what
  the reference MATLAB implementation optimizes: its score divides by a
  `2 * lambda` factor that is constant within a call and so cannot
  influence the arg max, leaving `lambda` inert. Provided for
  reproducing published results.

- `"density"`:

  \\w/m\\. Degenerates toward the single densest pair of nodes and is
  only useful with a large `min_size`.

## References

Chen S, Zhang Y, Wu Q, Bi C, Kochunov P, Hong LE (2023). Identifying
covariate-related subnetworks for whole-brain connectome analysis.
*Biostatistics*, kxad007.
[doi:10.1093/biostatistics/kxad007](https://doi.org/10.1093/biostatistics/kxad007)

Wu Q, Huang X, Culbreth AJ, Waltz JA, Hong LE, Chen S (2022). Extracting
brain disease-related connectome subgraphs by adaptive dense subgraph
discovery. *Biometrics* 78(4), 1566-1578.
[doi:10.1111/biom.13537](https://doi.org/10.1111/biom.13537)

## See also

[`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md) for
extraction plus inference,
[`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md)
for choosing `lambda` and `threshold`.

## Examples

``` r
sim <- simulate_fc(n = 100, n_nodes = 60, cluster_size = 12, f2 = 0.15,
                   seed = 1)
part <- subnet_extract(sim$W, threshold = quantile(vech(sim$W), 0.95))
part
#> Greedy-peeling partition
#>   nodes           : 60
#>   objective       : sicers (lambda = 0.60)
#>   threshold       : 2.537  (5.0% of edges retained)
#>   subnetworks     : 8
#> 
#>  subnet size edges density
#>       1   12    52   0.788
#>       2    3     2   0.667
#>       3    3     3   1.000
#>       4    5     3   0.300
#>       5    4     2   0.333
#>       6    4     2   0.333
#>       7    4     2   0.333
#>       8    4     2   0.333
#> 
#>   background      : 21 nodes
```
