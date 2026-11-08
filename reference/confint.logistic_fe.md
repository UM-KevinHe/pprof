# Confidence intervals for the provider effects of a logistic fixed-effect fit

Intervals for \\\gamma_i\\ obtained by inverting a test: the exact
Poisson-binomial test, the score test with \\\hat\beta\\ held fixed, or
a Wald test with the naive or patient-clustered standard error.
Intervals for standardized measures are produced by
[`standardize()`](https://um-kevinhe.github.io/pprof/reference/standardize.md),
which maps these limits through the standardization.

## Usage

``` r
# S3 method for class 'logistic_fe'
confint(
  object,
  parm = NULL,
  level = 0.95,
  type = c("exact", "score", "wald", "robust_wald"),
  wald_variance = c("case_mix", "full"),
  alternative = c("two.sided", "greater", "less"),
  ...
)
```

## Arguments

- object:

  a model fitted by
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md).

- parm:

  optional subset of provider identifiers.

- level:

  confidence level. Default 0.95.

- type:

  `"exact"` (default), `"score"`, `"wald"` or `"robust_wald"`.

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

## Details

Providers with no events have no finite lower limit and providers with
all events no finite upper limit; their finite limit is one-sided at
level `1 - level`. The exact inversion works on the individual Bernoulli
trials and is slow for providers with very many trials.

## Examples

``` r
d <- sim_kecc(n_providers = 10, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_fe(pd)
confint(fit, type = "score")
#>    provider   estimate      lower       upper
#> 1      P001 -0.5016517 -1.0249137  0.02184901
#> 2      P002 -0.9618830 -1.2923892 -0.63132908
#> 3      P003 -0.1930757 -0.7812240  0.39521893
#> 4      P004 -1.4712752 -1.8799614 -1.06247731
#> 5      P005 -0.9757544 -1.3901379 -0.56124741
#> 6      P006 -0.6629494 -0.9685020 -0.35736727
#> 7      P007 -0.4959592 -0.9180087 -0.07386339
#> 8      P008 -2.0741234 -2.6519146 -1.49608335
#> 9      P009 -0.6788108 -1.0268220 -0.33076124
#> 10     P010 -0.6808666 -1.0730397 -0.28860917
```
