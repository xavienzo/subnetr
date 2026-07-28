# Changelog

## subnetr 0.1.0

First release. An R reimplementation of the `Subnet` MATLAB tool for
detecting predictor-associated functional connectivity subnetworks, with
simulation-based power analysis.

### Analysis

- [`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md)
  runs the full pipeline: screen, tune, extract, test.
- [`subnet_extract()`](https://xavienzo.github.io/subnetr/reference/subnet_extract.md)
  exposes greedy peeling on its own.
- [`edge_stats()`](https://xavienzo.github.io/subnetr/reference/edge_stats.md)
  computes covariate-adjusted edge-wise association statistics from
  subject-level connectivity, vectorized over edges via the
  Frisch-Waugh-Lovell theorem and evaluated on the log p-value scale so
  that extreme evidence does not underflow.
- [`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md)
  selects the screening threshold and the objective parameter.
- Reordered-matrix and power-curve plots, plus
  [`membership()`](https://xavienzo.github.io/subnetr/reference/membership.md),
  [`dice()`](https://xavienzo.github.io/subnetr/reference/dice.md),
  [`vech()`](https://xavienzo.github.io/subnetr/reference/vech.md) /
  [`unvech()`](https://xavienzo.github.io/subnetr/reference/vech.md).

### Power analysis

- [`power_curve()`](https://xavienzo.github.io/subnetr/reference/power_curve.md)
  estimates detection and recovery power over a grid of sample sizes;
  [`required_n()`](https://xavienzo.github.io/subnetr/reference/required_n.md)
  reads a target off the curve.
- [`simulate_fc()`](https://xavienzo.github.io/subnetr/reference/simulate_fc.md)
  plants subnetworks with specified effect size and contamination rates.
  Its default `method = "fast"` draws edge t-statistics from their exact
  non-central sampling distribution instead of materializing
  subject-level data, which is what makes large replicate counts
  practical.

### Differences from the reference MATLAB implementation

- **`lambda` now has an effect.** `greedy_peeling_v2.m` divides its
  score by `2 * lambda`, a constant within a call, so `lambda` could not
  influence the arg max and the `lambda` search in `param_tuning.m` was
  a no-op. The generalized objective of Chen et al. (2023) is restored
  as the default; the original behaviour remains available as
  `objective = "avg_degree"` and is checked for exact agreement in the
  test suite.
- **Tuning no longer inflates the type-I error.** Selecting the
  threshold from the data and then testing against a null built only at
  the selected value roughly doubles the error rate.
  [`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md)
  defaults to `null = "retune"`, which makes every permutation run the
  same parameter search, restoring calibration at no extra asymptotic
  cost.
- **P-values use the add-one estimator**
  `(1 + #{T_null <= T_obs}) / (M + 1)`, so they are never exactly zero.
- **Extraction stops when the remainder carries no supra-threshold
  weight** rather than continuing until `N - 1` nodes are assigned.

### Performance

- The screened graph is stored sparsely and peeled with a lazily-updated
  min-heap, so a peeling pass costs `O(E log E)` rather than `O(N^2)`.
- A permutation places only the surviving edges, `O(E)` rather than
  `O(N^2)`, using the fact that permuting a screened weight vector is
  equivalent in distribution to choosing `E` positions at random and
  dealing the surviving weights into them.
- End to end this is roughly 8 to 14 times faster than a vectorized R
  port of the reference algorithm. See `benchmarks/benchmark.R`.
- Permutations seed their own RNG streams from their absolute index, so
  results depend on `seed` alone and never on `n_cores` or on how work
  was chunked.
