# Plot the reordered connectivity matrix

Draws the edge-weight matrix with nodes reordered so that the extracted
subnetworks sit in consecutive blocks along the diagonal, and outlines
the significant blocks. This is the plot to look at when judging whether
a detection is a compact subnetwork or a diffuse smear.

## Usage

``` r
# S3 method for class 'subnet'
plot(
  x,
  what = c("reordered", "observed"),
  significant_only = TRUE,
  col = grDevices::hcl.colors(64, "Inferno"),
  main = NULL,
  ...
)
```

## Arguments

- x:

  A `subnet` object.

- what:

  `"reordered"` (default) or `"observed"` for the matrix in its original
  node order.

- significant_only:

  Outline only the significant subnetworks.

- col:

  Colour palette.

- main:

  Plot title.

- ...:

  Passed to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html).

## Value

`x`, invisibly. Called for the plot.

## Examples

``` r
sim <- simulate_fc(n = 120, n_nodes = 60, cluster_size = 12, f2 = 0.2,
                   seed = 11)
fit <- subnet(sim$W, n_perm = 99, seed = 1)
plot(fit)

```
