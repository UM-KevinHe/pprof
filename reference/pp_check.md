# Data checks for a provider model

Missing values, (near-)zero variance of covariates, pairwise correlation
among covariates and variance inflation factors.
[`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
runs these checks when `check = TRUE`; the function is exported for
exploratory use.

## Usage

``` r
pp_check(
  data,
  y,
  x,
  provider,
  n = NULL,
  cor_threshold = 0.9,
  vif_threshold = 10,
  message = TRUE
)
```

## Arguments

- data:

  a data frame.

- y:

  name of the outcome column: 0/1 for a binary outcome, the number of
  events out of `n` trials for a binomial outcome, or any numeric value
  for a continuous outcome. `NULL` for a survival outcome, which is
  given by `time` and `event`.

- x:

  names of the covariate columns. Covariates must be numeric; code
  factors as 0/1 columns first.

- provider:

  name of the provider identifier column.

- n:

  name of the column with the number of trials per row (binomial
  outcome), or `NULL`.

- cor_threshold:

  warn about pairs of covariates whose absolute correlation exceeds this
  value. Default 0.9.

- vif_threshold:

  warn about covariates whose variance inflation factor is at least this
  value. Default 10.

- message:

  print a short summary. Default `TRUE`.

## Value

Invisibly, a list: `missing` (missing values per column), `nzv`
(near-zero-variance metrics per covariate), `cor` (correlation matrix),
`high_cor` (pairs above the threshold) and `vif`.

## Details

Missing values in any named column and covariates with zero variance are
errors. Near-zero variance (frequency ratio of the two most common
values above 95/5 and at most 10% distinct values), high correlation and
large VIF are warnings. The VIF of a covariate is \\1/(1 - R^2)\\ from
regressing it on the other covariates.

## Examples

``` r
d <- sim_kecc(n_providers = 10, seed = 2)
chk <- pp_check(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider", message = FALSE)
chk$vif
#>       x1       x2       x3 
#> 1.009384 1.002768 1.006736 
```
