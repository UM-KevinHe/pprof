# Model-based variance of a linear fixed-effect fit

Model-based variance of a linear fixed-effect fit

## Usage

``` r
# S3 method for class 'linear_fe'
vcov(object, ...)
```

## Arguments

- object:

  a model fitted by
  [`fit_linear_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_linear_fe.md).

- ...:

  unused.

## Value

The covariance matrix of `beta`, with the `"full"` variances of the
provider effects (see
[`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md))
as the attribute `"var_gamma"`.
