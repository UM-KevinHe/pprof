# Group providers by size

Group providers by size

## Usage

``` r
assign_groups(size, n_groups = 4, grouping = c("rank", "quantile", "cut"))
```

## Arguments

- size:

  one size per provider.

- n_groups:

  number of groups.

- grouping:

  `"rank"`, `"quantile"` or `"cut"`, see
  [`null_model()`](https://um-kevinhe.github.io/pprof/reference/null_models.md).

## Value

Integer group labels `1, ..., n_groups`.

## Examples

``` r
assign_groups(c(5, 50, 12, 30, 8, 100), n_groups = 3)
#> [1] 1 3 2 2 1 3
```
