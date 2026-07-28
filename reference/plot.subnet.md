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
  weight_label = NULL,
  legend = TRUE,
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

- weight_label:

  Label for the colour scale, as a string or an
  [expression](https://rdrr.io/r/base/expression.html). Defaults to the
  `"statistic"` attribute of `W`, and to `"edge weight"` when the matrix
  carries none.

- legend:

  Draw the colour scale.

- ...:

  Passed to
  [`graphics::image()`](https://rdrr.io/r/graphics/image.html).

## Value

`x`, invisibly. Called for the plot.

## Details

The colour scale is labelled with the quantity being displayed. When `W`
came from
[`edge_stats()`](https://xavienzo.github.io/subnetr/reference/edge_stats.md)
or
[`simulate_fc()`](https://xavienzo.github.io/subnetr/reference/simulate_fc.md)
that label is carried along with the matrix, so the figure states
whether it is showing \\-\log\_{10}\\ p, a t statistic, or something
else. Supply `weight_label` for matrices built by other means.

## Examples

``` r
sim <- simulate_fc(n = 120, n_nodes = 60, cluster_size = 12, f2 = 0.2,
                   seed = 11)
fit <- subnet(sim$W, n_perm = 99, seed = 1)
plot(fit)


# A matrix of t statistics from a real study labels itself accordingly
plot(fit, weight_label = "t statistic")

```
