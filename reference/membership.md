# Subnetwork membership for every node

Subnetwork membership for every node

## Usage

``` r
membership(x, significant_only = TRUE)
```

## Arguments

- x:

  A `subnet` object.

- significant_only:

  Label only the subnetworks that reached significance, leaving the rest
  as background.

## Value

An integer vector of length `n_nodes`. Entry `i` is the index of the
subnetwork containing node `i`, or `0` for background nodes.

## Examples

``` r
sim <- simulate_fc(n = 100, n_nodes = 60, cluster_size = 12, f2 = 0.2,
                   seed = 5)
table(membership(subnet(sim$W, n_perm = 99, seed = 1)))
#> 
#>  0  1 
#> 48 12 
```
