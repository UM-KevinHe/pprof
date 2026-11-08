# Standardized differences of a linear fixed-effect fit

Indirect or direct standardization of every provider's mean outcome, as
a difference from the expected mean.

## Usage

``` r
# S3 method for class 'linear_fe'
standardize(
  fit,
  type = c("indirect", "direct"),
  null = "median",
  providers = NULL,
  gamma_se = NULL,
  level = NULL,
  ci = c("wald", "robust_wald"),
  wald_variance = c("case_mix", "full"),
  alternative = c("two.sided", "greater", "less"),
  ...
)
```

## Arguments

- fit:

  a model fitted by
  [`fit_linear_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_linear_fe.md).

- type:

  `"indirect"` (default) or `"direct"`.

- null:

  for indirect standardization, the norm of the provider effect:
  `"median"` (default), `"mean"` (weighted by rows) or a number.

- providers:

  optional subset of provider identifiers.

- gamma_se:

  optional standard errors of the provider effects, named by provider
  (for instance the `se` column of a
  [`test()`](https://um-kevinhe.github.io/pprof/reference/test.md)
  result); with them `se_difference` is returned.

- level:

  if given, confidence limits of the difference are added, mapped from
  Wald limits of \\\gamma_i\\.

- ci:

  the standard error of the limits: `"wald"` (default; model-based, with
  `wald_variance`) or `"robust_wald"`, see
  [`confint.linear_fe()`](https://um-kevinhe.github.io/pprof/reference/confint.linear_fe.md).

- wald_variance:

  `"case_mix"` (default) or `"full"`, see
  [`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md).

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

- ...:

  unused.

## Value

A data frame with one row per provider: `provider`, `observed`,
`expected`, `difference`, and when requested `se_difference`, `lower`,
`upper`. Attributes `type`, `null_value`, `level`.

## Details

With \\\mu\_{ij}(\gamma) = \gamma + Z\_{ij}^\top\hat\beta\\:

- indirect: `observed` \\O_i = \sum\_{j \in i} Y\_{ij}\\, `expected`
  \\E_i = \sum\_{j \in i} \mu\_{ij}(\gamma_0)\\, difference \\(O_i -
  E_i) / n_i = \hat\gamma_i - \gamma_0\\;

- direct: `expected` \\E_i = \sum_k \mu_k(\hat\gamma_i)\\ over every
  fitted row, `observed` \\O = \sum_k Y_k\\, difference \\(E_i - O) /
  n\\, the deviation of \\\hat\gamma_i\\ from the row-weighted mean of
  the effects.

Both differences are linear in \\\hat\gamma_i\\, so their standard error
with \\\hat\beta\\ held fixed is \\se(\hat\gamma_i)\\ and the limits are
those of \\\gamma_i\\ shifted.

## Examples

``` r
d <- sim_kecc(n_providers = 20, type = "continuous", seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_linear_fe(pd)
head(standardize(fit, level = 0.95))
#>   provider    observed   expected difference       lower      upper
#> 1     P001  -42.582888 -10.896083 -0.4951063 -0.76938850 -0.2208242
#> 2     P002 -125.896359 -67.439053 -0.2967376 -0.45307196 -0.1404032
#> 3     P003    8.964691  -8.723502  0.3930710  0.06597039  0.7201715
#> 4     P004   34.058194 -43.422138  0.4724410  0.30109834  0.6437838
#> 5     P005   41.090638 -17.888160  0.4834328  0.28477404  0.6820915
#> 6     P006 -101.501224 -39.079480 -0.3001045 -0.45224891 -0.1479602
```
