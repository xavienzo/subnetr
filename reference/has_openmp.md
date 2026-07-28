# Report whether the compiled code has OpenMP support

The permutation loop is parallelized with OpenMP where the toolchain
provides it. When it does not – most notably the default Apple clang on
macOS – `subnet_test()` falls back to forked R workers, so the `n_cores`
argument still delivers a speedup.

## Usage

``` r
has_openmp()
```

## Value

A single logical value.

## Examples

``` r
has_openmp()
#> [1] TRUE
```
