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
#> 2      2    8     9 0.3214286   -6.748776    1.00       FALSE
#> 3      3    5     4 0.4000000   -4.358685    1.00       FALSE
#> 4      4   14    14 0.1538462   -2.683971    1.00       FALSE
#> 5      5    5     3 0.3000000   -2.656538    1.00       FALSE
#> 6      6    4     2 0.3333333   -2.169235    1.00       FALSE
#> 7      7    3     2 0.6666667   -3.575551    1.00       FALSE
```
