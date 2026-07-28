# Plot a power curve

Plot a power curve

## Usage

``` r
# S3 method for class 'subnet_power'
plot(x, which = c("both", "any", "recovery"), ci = TRUE, target = 0.8, ...)
```

## Arguments

- x:

  A `subnet_power` object.

- which:

  Curves to draw: `"both"` (default), `"any"` or `"recovery"`.

- ci:

  Draw pointwise Monte Carlo 95% intervals.

- target:

  Optional horizontal reference line, e.g. `0.8`.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly. Called for the plot.

## Examples

``` r
# \donttest{
pw <- power_curve(n = c(40, 80, 120), n_nodes = 60, cluster_size = 12,
                  f2 = 0.2, n_sim = 20, n_perm = 99, seed = 1)
#> n = 40  (20 replicates)
#> n = 80  (20 replicates)
#> n = 120  (20 replicates)
plot(pw, target = 0.8)

# }
```
