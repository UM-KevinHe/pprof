# Covariance of the covariate effects of a Cox fit

Covariance of the covariate effects of a Cox fit

## Usage

``` r
# S3 method for class 'cox'
vcov(object, variance = c("naive", "robust"), stage = 1, ...)
```

## Arguments

- object:

  a model fitted by
  [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md).

- variance:

  `"naive"` (default; the inverse information) or `"robust"` (the
  patient-clustered sandwich, available when
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  has a patient column).

- stage:

  `1` (default) or `2`, the stage whose coefficients are wanted.

- ...:

  unused.

## Value

A covariance matrix.
