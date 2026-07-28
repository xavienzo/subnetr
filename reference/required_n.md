# Sample size required to reach a target power

Interpolates a
[`power_curve()`](https://xavienzo.github.io/subnetr/reference/power_curve.md)
to report the smallest sample size reaching a target power.

## Usage

``` r
required_n(object, target = 0.8, which = c("recovery", "any"))
```

## Arguments

- object:

  A `subnet_power` object.

- target:

  Target power, e.g. `0.8`.

- which:

  Which power curve to read: `"recovery"` (default, the probability of
  recovering the planted subnetworks) or `"any"` (the probability of any
  detection).

## Value

A single number: the interpolated sample size, `NA` with a warning if
the curve never reaches `target` over the range simulated.

## Details

Linear interpolation between the two bracketing sample sizes. Extend the
grid rather than trusting an extrapolation if the target lies outside
it.

## Examples

``` r
# \donttest{
pw <- power_curve(n = c(40, 80, 120), n_nodes = 60, cluster_size = 12,
                  f2 = 0.2, n_sim = 20, n_perm = 99, seed = 1)
#> n = 40  (20 replicates)
#> n = 80  (20 replicates)
#> n = 120  (20 replicates)
required_n(pw, target = 0.8)
#> [1] 40
# }
```
