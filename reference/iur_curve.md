# Reliability at representative provider sizes

Evaluates the fitted reliability curve \\s^2_b / (s^2_b + s2\\within /
n)\\ at the smallest size, the mean size of each quantile group of sizes
and the largest size.

## Usage

``` r
iur_curve(x, sizes = NULL, n_quantiles = 10)
```

## Arguments

- x:

  a `pp_iur` object from
  [`iur()`](https://um-kevinhe.github.io/pprof/reference/iur.md).

- sizes:

  sizes to group; default the providers' sizes.

- n_quantiles:

  number of quantile groups. Default 10.

## Value

A data frame with `size_group` (`"min"`, `"quantile 1"`, ..., `"max"`),
`size` and `iur`.

## Examples

``` r
r <- iur("direct", estimate = c(0.2, 0.3, 0.25, 0.4), se = c(0.05, 0.04, 0.06, 0.03),
                size = c(40, 60, 30, 90))
iur_curve(r, n_quantiles = 2)
#>   size_group size       iur
#> 1        min   30 0.6914213
#> 2 quantile 1   35 0.7233066
#> 3 quantile 2   75 0.8485230
#> 4        max   90 0.8704998
```
