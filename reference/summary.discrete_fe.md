# Tests of the covariates of a discrete-time fixed-effect fit

Wald tests of the covariate coefficients with the model-based variance
or the sandwich clustered by patient (by subject without a patient
column), the provider effects profiled out and the period effects
estimated jointly.

## Usage

``` r
# S3 method for class 'discrete_fe'
summary(
  object,
  test = c("wald", "robust_wald"),
  covariates = NULL,
  level = 0.95,
  null = 0,
  ...
)
```

## Arguments

- object:

  a model fitted by
  [`fit_discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_discrete_fe.md).

- test:

  `"wald"` (default) or `"robust_wald"`.

- covariates:

  optional subset of covariate names.

- level:

  confidence level. Default 0.95.

- null:

  the null value of the coefficients. Default 0.

- ...:

  unused.

## Value

A data frame of class `pp_summary` with `covariate`, `estimate`, `se`,
`statistic`, `p`, `lower`, `upper`.
