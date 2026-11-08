# Standardized measures of a random-effect fit

Indirect or direct standardization of every provider's outcome from a
random-effect model: ratios and rates for the logistic models,
differences for the linear models.

## Usage

``` r
# S3 method for class 'logistic_re'
standardize(
  fit,
  type = c("indirect", "direct"),
  numerator = c("observed", "predicted"),
  null = "median",
  n_pop = NULL,
  providers = NULL,
  alpha_se = NULL,
  level = NULL,
  alternative = c("two.sided", "greater", "less"),
  threads = 1,
  ...
)

# S3 method for class 'linear_re'
standardize(
  fit,
  type = c("indirect", "direct"),
  numerator = c("observed", "predicted"),
  null = "median",
  providers = NULL,
  alpha_se = NULL,
  level = NULL,
  alternative = c("two.sided", "greater", "less"),
  ...
)
```

## Arguments

- fit:

  a model from
  [`fit_logistic_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md),
  [`fit_logistic_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md),
  [`fit_linear_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md)
  or
  [`fit_linear_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md).

- type:

  `"indirect"` (default) or `"direct"`.

- numerator:

  for indirect standardization, the provider's `"observed"` outcome
  (default) or its `"predicted"` outcome under the fitted model with the
  provider's random effect, \\\sum\_{j \in i} N\_{ij}\\ \hat p\_{ij}\\
  (a shrunken numerator, as in predicted-to-expected ratios).

- null:

  for indirect standardization, the norm of the random effect:
  `"median"` (default) or `"mean"` (weighted by rows) of the predicted
  effects, or a number (0 is the population intercept).

- n_pop:

  denominator of the rates (logistic models); defaults to the number of
  trials among the fitted rows.

- providers:

  optional subset of provider identifiers.

- alpha_se:

  optional standard errors of the random effects, named by provider;
  with them the direct rate's standard error (logistic) or the
  difference's standard error (linear) is returned.

- level:

  if given, confidence limits are added, mapped from Wald limits of
  \\\alpha_i\\ (see
  [`confint.logistic_re()`](https://um-kevinhe.github.io/pprof/reference/confint_mixed.md));
  for the logistic models only with `numerator = "predicted"`, whose
  limits are the numerator evaluated at the limits of \\\alpha_i\\.

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

- threads:

  OpenMP threads for direct standardization.

- ...:

  unused.

## Value

A data frame with one row per provider. Logistic models: `provider`,
`observed` (the numerator), `expected`, `ratio`, `rate`,
`variance_expected` and, when requested, `se_rate`, `lower_ratio`,
`upper_ratio`, `lower_rate`, `upper_rate`. Linear models: `provider`,
`observed`, `expected`, `difference` and, when requested,
`se_difference`, `lower`, `upper`. Attributes `type`, `numerator`,
`null_value`, `level` (and `n_pop`, `population_rate` for the logistic
models).

## Details

With \\\eta\_{ij} = \hat\mu + Z\_{ij}^\top\hat\beta\\ the fixed part of
the linear predictor and \\p\_{ij}(a) = \mathrm{expit}(a +
\eta\_{ij})\\:

- logistic, indirect: \\E_i = \sum\_{j \in i} N\_{ij}\\
  p\_{ij}(\alpha_0)\\, ratio = numerator / \\E_i\\;

- logistic, direct: \\E_i = \sum_k N_k\\ p_k(\hat\alpha_i)\\ over every
  fitted row, ratio \\E_i / \sum O\\;

- linear, indirect: \\E_i = \sum\_{j \in i} (\alpha_0 + \eta\_{ij})\\,
  difference = (numerator \\- E_i) / n_i\\;

- linear, direct: \\E_i = \sum_k (\hat\alpha_i + \eta_k)\\, difference
  \\(E_i - \sum O) / n\\.

Rates are ratio \\\times \sum O / n\_{pop}\\, truncated to \[0, 1\].

## Examples

``` r
# \donttest{
d <- sim_kecc(n_providers = 20, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_re(pd)
head(standardize(fit, numerator = "predicted", level = 0.95))
#>   provider observed expected     ratio      rate variance_expected lower_ratio
#> 1     P001 12.56566 16.16176 0.7774934 0.1956068         11.308641   0.5063884
#> 2     P002 32.73440 43.73419 0.7484854 0.1883088         31.299285   0.5598575
#> 3     P003 18.18704 11.23172 1.6192565 0.4073829          7.728969   1.1844639
#> 4     P004 60.53738 38.90918 1.5558635 0.3914341         27.248355   1.2866897
#> 5     P005 33.76193 31.07871 1.0863364 0.2733075         21.580101   0.8303898
#> 6     P006 43.22211 51.87075 0.8332656 0.2096383         35.856866   0.6490632
#>   upper_ratio lower_rate upper_rate
#> 1   1.1385243  0.1274005  0.2864373
#> 2   0.9816633  0.1408525  0.2469731
#> 3   2.0921356  0.2979950  0.5263529
#> 4   1.8439836  0.3237137  0.4639212
#> 5   1.3830801  0.2089148  0.3479642
#> 6   1.0509853  0.1632955  0.2644136
# }
```
