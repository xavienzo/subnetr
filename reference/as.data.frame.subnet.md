# Tabulate detected subnetworks

Tabulate detected subnetworks

## Usage

``` r
# S3 method for class 'subnet'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  A `subnet` object.

- row.names, optional:

  Ignored, present for method consistency.

- ...:

  Ignored.

## Value

A data frame with one row per extracted subnetwork: its index, node
count, supra-threshold edge count, within-block density, log binomial
statistic, permutation p-value and significance flag.

## Examples

``` r
sim <- simulate_fc(n = 100, n_nodes = 60, cluster_size = 12, f2 = 0.2,
                   seed = 5)
as.data.frame(subnet(sim$W, n_perm = 99, seed = 1))
#>   subnet size edges   density   statistic p_value significant
#> 1      1   12    59 0.8939394 -116.103772    0.01        TRUE
#> 2      2   26    39 0.1200000   -2.005758    1.00       FALSE
#> 3      3    8     8 0.2857143   -5.308157    1.00       FALSE
#> 4      4    7     4 0.1904762   -1.884103    1.00       FALSE
```
