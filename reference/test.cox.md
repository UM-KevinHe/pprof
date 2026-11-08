# Tests of the standardized ratios of a Cox fit

Tests every provider's indirect standardized ratio \\O_j / E_j\\ (see
[`standardize.cox()`](https://um-kevinhe.github.io/pprof/reference/standardize.cox.md))
against 1, treating the observed count \\O_j\\ as Poisson with mean
\\E_j\\, and flags the providers with more or fewer events than
expected. These are the tests of the survival measures (SMR, SHR, STrR,
SWR).

## Usage

``` r
# S3 method for class 'cox'
test(
  fit,
  statistic = c("midp", "exact", "wald"),
  dispersion = c("none", "deviance"),
  providers = NULL,
  level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  min_expected = 0,
  null_model = NULL,
  ...
)
```

## Arguments

- fit:

  a model fitted by
  [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md).

- statistic:

  `"midp"` (default; the two-sided mid-p Poisson test turned into a z
  statistic, the SMR convention), `"exact"` (the two-sided Poisson test
  with exact or Byar limits, the SWR convention; theoretical null only)
  or `"wald"` (\\\log(O_j / E_j)\\ over its standard error \\\sqrt{\phi
  / O_j}\\, the SHR convention).

- dispersion:

  for `statistic = "wald"`: `"none"` (default; \\\phi = 1\\) or
  `"deviance"` (\\\phi\\ = the Poisson deviance of the patient-level
  counts about \\(O_j / E_j) e_k\\ over its degrees of freedom, as
  `PROC GENMOD`'s `DSCALE`; the units are the patients when
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  has a patient column, the rows otherwise).

- providers:

  optional subset of provider identifiers to report.

- level:

  confidence level; the significance level is `1 - level`.

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

- min_expected:

  providers with fewer than `min_expected` expected events are dropped
  from the result (the measures use 3, so that an empirical null fitted
  on the result does not see them). Default 0.

- null_model:

  optional null model: a `pp_null` object or a function of the result
  returning one, for instance
  `function(t) null_model("empirical", t, size = t$person_time, grouping = "quantile")`
  (the survival measures group providers by quartiles of person-time);
  the result is then passed through
  [`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md).
  Not available for the exact test.

- ...:

  unused.

## Value

A data frame of class `pp_test` with one row per tested provider:
`provider`, `estimate` (the ratio \\O_j / E_j\\), `se` (Wald only; the
standard error of the log ratio), `null_value` (1), `statistic`, `p`,
`flag` (`1` more events than expected, `-1` fewer, `0` not different,
`NA` not testable), `lower`, `upper` (limits for the ratio), `observed`,
`expected`, `person_time`. Attributes `test` (the settings, including
`phi`) and `sizes`.

## Details

The mid-p statistic is \\z = \pm\Phi^{-1}(p)\\ with \\p = \max(10^{-6},
\min(P\_{low}, P\_{high}) / 2)\\, \\P\_{low} = 2 F(O_j) - f(O_j)\\,
\\P\_{high} = 2 \[1 - F(O_j - 1)\] - f(O_j)\\, \\F\\ and \\f\\ the
Poisson distribution and mass functions with mean \\E_j\\, the sign
negative when the lower tail is the smaller. Its confidence limits
invert the (calibrated) test: the Poisson means at which the provider's
count is just significant, divided by \\E_j\\, so that a limit excludes
1 exactly when the provider is flagged.

The exact test's p-value is \\2 P(X \ge O_j)\\ for a ratio above 1 and
\\2 P(X \le O_j)\\ otherwise, capped at 0.999; its limits are the exact
chi-square limits when \\E_j \< 100\\ and Byar's approximation
otherwise.

The Wald limits are \\\exp\[\log(O_j / E_j) \pm c\\ se\]\\; a provider
with no events has no Wald statistic.

## References

He K, Schaubel DE (2014). Standardized mortality ratio for evaluating
center-specific mortality: assessment and alternative. *Statistics in
Biosciences*, 6(2), 1-19.

Breslow NE, Day NE (1987). *Statistical Methods in Cancer Research,
Volume II: The Design and Analysis of Cohort Studies*. IARC.

## Examples

``` r
s <- sim_kecc(n_providers = 30, type = "survival", seed = 1)
ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", time = "stop", start = "start",
              event = "event", cutoff = 0, message = FALSE)
fit <- fit_cox(ps)
test(fit)
#> Provider tests: midp statistic, poisson variance; null ratio = 1; two.sided, level 0.95
#>   30 providers: 3 above, 4 below, 23 not different, 0 not tested
#>  provider estimate se null_value statistic        p flag  lower upper observed
#>      P001   4.5464 NA          1    4.7534 0.000002    1 2.6913 7.226       16
#>      P002   1.2270 NA          1    1.0568 0.290581    0 0.8251 1.761       27
#>      P003   1.1627 NA          1    0.4457 0.655800    0 0.5085 2.300        7
#>      P004   0.7691 NA          1   -1.0881 0.276549    0 0.4630 1.206       17
#>      P005   1.4076 NA          1    1.4777 0.139493    0 0.8840 2.135       20
#>      P006   1.4212 NA          1    1.8140 0.069676    0 0.9699 2.015       29
#>  expected person_time
#>     3.519       6.915
#>    22.005      48.263
#>     6.020      13.663
#>    22.104      45.257
#>    14.209      31.364
#>    20.405      45.432
#>   ... 24 more rows
test(fit, statistic = "wald", dispersion = "deviance")
#> Provider tests: wald statistic, poisson variance; null ratio = 1; two.sided, level 0.95
#>   30 providers: 2 above, 3 below, 25 not different, 0 not tested
#>  provider estimate     se null_value statistic         p flag  lower upper
#>      P001   4.5464 0.2848          1    5.3169 1.056e-07    1 2.6015 7.945
#>      P002   1.2270 0.2193          1    0.9330 3.508e-01    0 0.7984 1.886
#>      P003   1.1627 0.4306          1    0.3501 7.263e-01    0 0.5000 2.704
#>      P004   0.7691 0.2763          1   -0.9502 3.420e-01    0 0.4475 1.322
#>      P005   1.4076 0.2547          1    1.3420 1.796e-01    0 0.8544 2.319
#>      P006   1.4212 0.2116          1    1.6615 9.662e-02    0 0.9388 2.151
#>  observed expected person_time
#>        16    3.519       6.915
#>        27   22.005      48.263
#>         7    6.020      13.663
#>        17   22.104      45.257
#>        20   14.209      31.364
#>        29   20.405      45.432
#>   ... 24 more rows
```
