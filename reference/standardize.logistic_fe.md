# Standardized ratios and rates of a logistic fixed-effect fit

Indirect or direct standardization of every provider's outcome, with the
observed and expected counts, the delta-method numerator of the direct
standard error, and optional confidence intervals mapped from the limits
of the provider effects.

## Usage

``` r
# S3 method for class 'logistic_fe'
standardize(
  fit,
  type = c("indirect", "direct"),
  null = "median",
  n_pop = NULL,
  providers = NULL,
  gamma_se = NULL,
  level = NULL,
  ci = c("exact", "score", "wald", "robust_wald"),
  alternative = c("two.sided", "greater", "less"),
  threads = 1,
  ...
)
```

## Arguments

- fit:

  a model fitted by
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md).

- type:

  `"indirect"` (default) or `"direct"`.

- null:

  for indirect standardization, the norm of the provider effect:
  `"median"` (default), `"mean"` (weighted by rows) or a number.

- n_pop:

  denominator of the rates: the number of trials in the standard
  population. Defaults to the number of trials among the fitted rows.

- providers:

  optional subset of provider identifiers.

- gamma_se:

  optional standard errors of the provider effects, named by provider
  (for instance the `se` column of a Wald
  [`test()`](https://um-kevinhe.github.io/pprof/reference/test.md)
  result); with them the direct rate's standard error
  \\\mathrm{numer\\se}/n\_{pop} \times se(\hat\gamma_i)\\ is returned.

- level:

  if given, confidence limits of the ratio and the rate are added,
  mapped from the limits of \\\gamma_i\\ obtained with `ci`.

- ci:

  the test inverted for the limits, see
  [`confint.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/confint.logistic_fe.md).

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

- threads:

  OpenMP threads for direct standardization.

- ...:

  unused.

## Value

A data frame with one row per provider: `provider`, `observed`,
`expected`, `ratio`, `rate`, `variance_expected` (indirect: \\\sum N
p_0(1 - p_0)\\; direct: `numer_se`), and when requested `se_rate`,
`lower_ratio`, `upper_ratio`, `lower_rate`, `upper_rate`. Attributes
`type`, `null_value`, `n_pop` and `population_rate`.

## Details

With \\p\_{ij}(\gamma) = \mathrm{expit}(\gamma + Z\_{ij}^\top
\hat\beta)\\ and \\N\_{ij}\\ trials per row:

- indirect: \\E_i = \sum\_{j \in i} N\_{ij}\\ p\_{ij}(\gamma_0)\\, ratio
  \\O_i / E_i\\;

- direct: \\E_i = \sum_k N_k\\ p_k(\hat\gamma_i)\\ over every fitted
  row, ratio \\E_i / \sum O\\.

In both cases rate = ratio \\\times \sum O / n\_{pop}\\, truncated to
\[0, 1\]. The direct `numer_se` is \\\sum_k N_k p_k(1 - p_k)\\ at
\\\hat\gamma_i\\; the direct rate's standard error with \\\hat\beta\\
held fixed is `numer_se / n_pop * se(gamma_i)`. For an arbitrary
standard population use
[`direct_expectation()`](https://um-kevinhe.github.io/pprof/reference/direct_expectation.md).

## References

He K (2019). Indirect and direct standardization for evaluating
transplant centers. *Journal of Hospital Administration*, 8(1), 9-14.

## Examples

``` r
d <- sim_kecc(n_providers = 20, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_fe(pd)
head(standardize(fit, type = "direct"))
#>   provider observed  expected     ratio      rate variance_expected
#> 1     P001       11  442.7491 0.6568978 0.1652666          342.9898
#> 2     P002       31  464.3512 0.6889483 0.1733300          355.4082
#> 3     P003       21 1219.3257 1.8090886 0.4551421          589.2703
#> 4     P004       63 1042.4641 1.5466826 0.3891243          567.2390
#> 5     P005       34  708.7526 1.0515617 0.2645586          472.3015
#> 6     P006       42  524.6180 0.7783650 0.1958260          388.2058
```
