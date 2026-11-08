# Model-based variance of a logistic fixed-effect fit

The inverse of the observed information for \\\beta\\ after profiling
out the provider effects. The model-based variances of the provider
effects are returned as the attribute `"var_gamma"`; see
[`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md)
for the choice of variance used by the Wald test.

## Usage

``` r
# S3 method for class 'logistic_fe'
vcov(object, ...)
```

## Arguments

- object:

  a model fitted by
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md).

- ...:

  unused.

## Value

A covariance matrix for `beta` with attribute `"var_gamma"`.
