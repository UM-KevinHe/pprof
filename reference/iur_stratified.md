# Reliability within strata of providers

Repeats the variance decomposition of
[`iur()`](https://um-kevinhe.github.io/pprof/reference/iur.md) within
strata of providers.

## Usage

``` r
iur_stratified(x, by = NULL, cuts = NULL)
```

## Arguments

- x:

  a `pp_iur` object from
  [`iur()`](https://um-kevinhe.github.io/pprof/reference/iur.md).

- by:

  one stratification value per provider; default the sizes.

- cuts:

  right-closed cut points of `by`; default its deciles.

## Value

A data frame with one row for all providers and one per stratum:
`stratum`, `n_providers`, `iur`.

## Examples

``` r
r <- iur("direct", estimate = c(0.2, 0.3, 0.25, 0.4, 0.1, 0.35),
                se = c(0.05, 0.04, 0.06, 0.03, 0.07, 0.03),
                size = c(40, 60, 30, 90, 20, 80))
iur_stratified(r, cuts = 50)
#>            stratum n_providers       iur
#> 1              all           6 0.8125807
#> 2 20 <= size <= 50           3 0.2803190
#> 3  50 < size <= 90           3 0.5517684
```
