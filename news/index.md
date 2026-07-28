# Changelog

## subnetr 0.1.0

First release. Detection of predictor-associated functional connectivity
subnetworks, with simulation-based power analysis for study planning.

### Analysis

- [`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md)
  runs the full pipeline: screen, tune, extract, test.
- [`subnet_extract()`](https://xavienzo.github.io/subnetr/reference/subnet_extract.md)
  exposes greedy peeling on its own. A candidate set of `n` nodes
  holding total weight `w` scores `w / n^(2 * lambda)`, a family that
  spans total weight at `lambda = 0`, average degree at `lambda = 0.5`,
  and edge density as `lambda` approaches 1.
- [`edge_stats()`](https://xavienzo.github.io/subnetr/reference/edge_stats.md)
  computes covariate-adjusted edge-wise association statistics from
  subject-level connectivity, vectorized over edges via the
  Frisch-Waugh-Lovell theorem and evaluated on the log p-value scale so
  that extreme evidence does not underflow.
- [`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md)
  selects the screening threshold and `lambda`, standardizing each
  candidate against its own permutation null so that settings producing
  graphs of different sparsity stay comparable.
- Inference is a Westfall-Young max-statistic permutation test, so
  reported p-values control the family-wise error rate across all
  extracted subnetworks with no further correction. They use the add-one
  estimator `(1 + #{T_null <= T_obs}) / (M + 1)` and are therefore never
  exactly zero.
- Choosing the threshold from the data and then testing against a null
  built only at the chosen value roughly doubles the type-I error.
  [`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md)
  defaults to `null = "retune"`, which has every permutation run the
  same parameter search, restoring calibration at no extra asymptotic
  cost.
- Reordered-matrix and power-curve plots, plus
  [`membership()`](https://xavienzo.github.io/subnetr/reference/membership.md),
  [`dice()`](https://xavienzo.github.io/subnetr/reference/dice.md),
  [`vech()`](https://xavienzo.github.io/subnetr/reference/vech.md) /
  [`unvech()`](https://xavienzo.github.io/subnetr/reference/vech.md).
  Matrix plots label their colour scale with the statistic being
  displayed, carried through from
  [`edge_stats()`](https://xavienzo.github.io/subnetr/reference/edge_stats.md).

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

### Performance

- The screened graph is stored sparsely and peeled with a lazily-updated
  min-heap, so a peeling pass costs `O(E log E)` rather than `O(N^2)` in
  the number of nodes.
- A permutation places only the surviving edges, `O(E)` rather than
  `O(N^2)`, using the fact that permuting a screened weight vector is
  equivalent in distribution to choosing `E` positions at random and
  dealing the surviving weights into them.
- A 200-permutation test on a 400-node connectome takes about 0.4 s
  single-core. See `benchmarks/benchmark.R`.
- Permutations seed their own RNG streams from their absolute index, so
  results depend on `seed` alone and never on `n_cores` or on how work
  was chunked.
