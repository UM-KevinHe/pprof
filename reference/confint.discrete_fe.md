# Confidence intervals for the provider effects of a discrete-time fit

Intervals for \\\gamma_i\\ by inverting the score test (naive or
clustered variance) or from a Wald test, as attached by
[`test.discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/test.discrete_fe.md).

## Usage

``` r
# S3 method for class 'discrete_fe'
confint(
  object,
  parm = NULL,
  level = 0.95,
  type = c("score", "robust_score", "wald", "robust_wald"),
  wald_variance = c("case_mix", "full"),
  alternative = c("two.sided", "greater", "less"),
  ...
)
```

## Arguments

- object:

  a model fitted by
  [`fit_discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_discrete_fe.md).

- parm:

  optional subset of provider identifiers.

- level:

  confidence level. Default 0.95.

- type:

  `"score"` (default), `"robust_score"`, `"wald"` or `"robust_wald"`.

- wald_variance:

  for `type = "wald"`: `"case_mix"` (default) or `"full"`, see
  [`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md).

- alternative:

  `"two.sided"` (default), `"greater"` (lower limit only) or `"less"`
  (upper limit only).

- ...:

  unused.

## Value

A data frame with `provider`, `estimate`, `lower`, `upper`.
