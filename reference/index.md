# Package index

## Analysis

Going from subject-level connectivity to significant subnetworks.

- [`subnet()`](https://xavienzo.github.io/subnetr/reference/subnet.md) :
  Detect predictor-associated subnetworks with a permutation test
- [`edge_stats()`](https://xavienzo.github.io/subnetr/reference/edge_stats.md)
  : Edge-wise association statistics for a connectivity study
- [`subnet_extract()`](https://xavienzo.github.io/subnetr/reference/subnet_extract.md)
  : Extract candidate subnetworks by greedy peeling
- [`tune_subnet()`](https://xavienzo.github.io/subnetr/reference/tune_subnet.md)
  : Tune the screening threshold and objective parameter

## Power analysis

Sizing a study before you run it.

- [`power_curve()`](https://xavienzo.github.io/subnetr/reference/power_curve.md)
  : Power analysis for subnetwork detection
- [`required_n()`](https://xavienzo.github.io/subnetr/reference/required_n.md)
  : Sample size required to reach a target power
- [`simulate_fc()`](https://xavienzo.github.io/subnetr/reference/simulate_fc.md)
  : Simulate a connectome study with planted predictor-associated
  subnetworks

## Inspecting results

- [`plot(`*`<subnet>`*`)`](https://xavienzo.github.io/subnetr/reference/plot.subnet.md)
  : Plot the reordered connectivity matrix
- [`plot(`*`<subnet_power>`*`)`](https://xavienzo.github.io/subnetr/reference/plot.subnet_power.md)
  : Plot a power curve
- [`as.data.frame(`*`<subnet>`*`)`](https://xavienzo.github.io/subnetr/reference/as.data.frame.subnet.md)
  : Tabulate detected subnetworks
- [`membership()`](https://xavienzo.github.io/subnetr/reference/membership.md)
  : Subnetwork membership for every node

## Utilities

- [`vech()`](https://xavienzo.github.io/subnetr/reference/vech.md)
  [`unvech()`](https://xavienzo.github.io/subnetr/reference/vech.md) :
  Convert between a symmetric matrix and its vectorized lower triangle
- [`dice()`](https://xavienzo.github.io/subnetr/reference/dice.md) :
  Sorensen-Dice overlap between two node sets
- [`has_openmp()`](https://xavienzo.github.io/subnetr/reference/has_openmp.md)
  : Report whether the compiled code has OpenMP support
