# subnetr: Detection and Power Analysis for Predictor-Associated Functional Connectivity Subnetworks

Identifies connected subnetworks of the functional connectome whose
edges are associated with a predictor of interest, and provides
simulation-based power analysis for study planning. Edge-level
association statistics are screened into a weighted graph, a generalized
densest-subgraph objective is optimized by greedy peeling to extract
candidate subnetworks, and family-wise error rate is controlled by a
max-statistic permutation test that accounts for the data-driven choice
of screening threshold. The peeling and permutation machinery is
implemented in C++ over a sparse representation of the screened graph,
which makes the repeated model fitting required for power analysis
practical. Methods follow Wu et al. (2022)
[doi:10.1111/biom.13537](https://doi.org/10.1111/biom.13537) and Chen et
al. (2023)
[doi:10.1093/biostatistics/kxad007](https://doi.org/10.1093/biostatistics/kxad007)
.

## See also

Useful links:

- <https://github.com/xavienzo/subnetr>

- Report bugs at <https://github.com/xavienzo/subnetr/issues>

## Author

**Maintainer**: Subnet Developers <yezhipan7@gmail.com>
