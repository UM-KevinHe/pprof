# Model-based variance of a discrete-time fixed-effect fit

The inverse of the observed information for the covariate coefficients
after profiling out the provider effects (the period effects are
estimated jointly with them). The model-based variances of the provider
effects are returned as the attribute `"var_gamma"`, those of the period
effects (periods 2 onwards) as `"var_alpha"`.

## Usage

``` r
# S3 method for class 'discrete_fe'
vcov(object, ...)
```

## Arguments

- object:

  a model fitted by
  [`fit_discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_discrete_fe.md).

- ...:

  unused.

## Value

A covariance matrix for `beta` with attributes `"var_gamma"` and
`"var_alpha"`.
