# Tests of the covariate coefficients of a Cox fit

Wald tests of the log hazard ratios of a stage of
[`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md).

## Usage

``` r
# S3 method for class 'cox'
summary(
  object,
  variance = c("naive", "robust"),
  covariates = NULL,
  level = 0.95,
  null = 0,
  stage = 1,
  ...
)
```

## Arguments

- object:

  a model fitted by
  [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md).

- variance:

  `"naive"` (default; the inverse information) or `"robust"` (the
  patient-clustered sandwich; needs a patient column).

- covariates:

  optional subset of covariate names.

- level:

  confidence level for the Wald intervals.

- null:

  the null value of each coefficient. Default 0.

- stage:

  `1` (default) or `2`.

- ...:

  unused.

## Value

A data frame with one row per covariate: `covariate`, `estimate` (the
log hazard ratio), `se`, `statistic` (z), `p`, `lower` and `upper`.

## Examples

``` r
s <- sim_kecc(n_providers = 20, type = "survival", seed = 1)
ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", time = "stop", start = "start",
              event = "event", cutoff = 0, message = FALSE)
summary(fit_cox(ps), variance = "robust")
#>   covariate   estimate         se statistic            p      lower     upper
#> 1        x1  0.5168166 0.05420841  9.533883 1.515082e-21  0.4105701 0.6230632
#> 2        x2 -0.1163487 0.11121845 -1.046128 2.955020e-01 -0.3343329 0.1016354
#> 3        x3  0.2665547 0.05777669  4.613534 3.958795e-06  0.1533145 0.3797949
```
