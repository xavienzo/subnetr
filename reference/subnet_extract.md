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
  objective = c("generalized", "density"),
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

  Size-penalty exponent in \[0, 1\]; see the section on choosing an
  objective. Ignored when `objective = "density"`.

- objective:

  One of `"generalized"` (default) or `"density"`.

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

A candidate node set \\S\\ of `n` nodes carrying total suprathreshold
edge weight `w` is scored by

- `"generalized"`:

  \\w / n^{2\lambda}\\, the default. This is the adaptive density
  function \\f(S; \lambda_W) = \|W(S)\| / \|S\|^{\lambda_W}\\ of Wu et
  al. (2022), and equivalently the \\\ell_0\\ graph norm shrinkage
  criterion \\\log\\U\\\_1 - \lambda_0 \log\\U\\\_0\\ of Chen et al.
  (2023), which rewards edge weight within a subnetwork while penalizing
  its size.

- `"density"`:

  \\w / \binom{n}{2}\\, pure edge density, ignoring `lambda`. With no
  size floor this degenerates to the single densest pair of nodes, so
  use it only with a substantial `min_size`.

## Parameterization of lambda

The `lambda` used here equals \\\lambda_0\\ of Chen et al. (2023) and
half of \\\lambda_W\\ of Wu et al. (2022), whose exponent runs over \[1,
2\]. Concretely, `lambda = 0.5` (\\\lambda_W = 1\\) gives the degree
density \\f_1\\, the objective of Charikar (2000), and `lambda = 1`
(\\\lambda_W = 2\\) gives the area density \\f_2\\. Larger values favour
smaller and denser subnetworks; `lambda = 0` places every node in one
subnetwork. Values between 0.5 and 0.7 are the usual working range, and
[`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md)
selects a value from the data when none is supplied.

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
#>   objective       : generalized (lambda = 0.60)
#>   threshold       : 2.537  (5.0% of edges retained)
#>   subnetworks     : 7
#> 
#>  subnet size edges density
#>       1   12    52   0.788
#>       2    6     7   0.467
#>       3    5     3   0.300
#>       4    4     2   0.333
#>       5    4     2   0.333
#>       6    4     2   0.333
#>       7    4     2   0.333
#> 
#>   background      : 21 nodes
```
