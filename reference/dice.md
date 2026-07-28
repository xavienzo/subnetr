# Sorensen-Dice overlap between two node sets

Sorensen-Dice overlap between two node sets

## Usage

``` r
dice(a, b)
```

## Arguments

- a, b:

  Integer vectors of node indices.

## Value

A number in \[0, 1\]; 1 means the two sets are identical and 0 that they
are disjoint. Returns 0 when either set is empty.

## Examples

``` r
dice(1:20, 11:30)
#> [1] 0.5
```
