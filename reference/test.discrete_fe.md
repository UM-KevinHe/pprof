# Tests of the provider effects of a discrete-time fixed-effect fit

Tests \\H_0: \gamma_i = \gamma_0\\ for every provider (or a subset) with
the score or the Wald statistic and a model-based or clustered variance,
and flags the providers above or below the null value.

## Usage

``` r
# S3 method for class 'discrete_fe'
test(
  fit,
  statistic = c("score", "wald"),
  variance = c("naive", "robust"),
  null = "median",
  providers = NULL,
  level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  wald_variance = c("case_mix", "full"),
  null_model = NULL,
  interval = TRUE,
  ...
)
```

## Arguments

- fit:

  a model fitted by
  [`fit_discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_discrete_fe.md).

- statistic:

  `"score"` (default) or `"wald"`.

- variance:

  `"naive"` (default; model-based) or `"robust"` (sandwich clustered by
  patient, or by subject when
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  has no patient column).

- null:

  the null value \\\gamma_0\\: `"median"` (default) or `"mean"`
  (weighted by person-periods) of the fitted effects, or a number.

- providers:

  optional subset of provider identifiers to report.

- level:

  confidence level; the significance level is `1 - level`.

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

- wald_variance:

  for `statistic = "wald"` with the naive variance: `"case_mix"`
  (default) or `"full"`, see
  [`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md).

- null_model:

  optional null model, as in
  [`test.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/test.logistic_fe.md);
  the result is then passed through
  [`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md).

- interval:

  for the score statistic, compute the confidence limits of \\\gamma_i\\
  by inverting the test (default `TRUE`).

- ...:

  unused.

## Value

A data frame of class `pp_test`, with the columns of
[`test.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/test.logistic_fe.md).

## Details

With \\H_j(g) = \sum\_{k \le t_j} \mathrm{expit}(g + \hat\alpha_k +
Z_j^\top\hat\beta)\\ the expected events of subject \\j\\ over its
periods at risk and \\W_j(g)\\ the sum of \\h(1 - h)\\:

- score: \\(O_i - \sum\_{j \in i} H_j(\gamma_0)) / \sqrt{V_i}\\ with the
  naive \\V_i = \sum\_{j \in i} W_j(\gamma_0)\\ or the clustered \\V_i =
  \sum_c (\sum\_{j \in c} \delta_j - H_j(\gamma_0))^2\\;

- Wald: \\(\hat\gamma_i - \gamma_0)/se\\, with the variance of
  [`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md)
  or the clustered standard error \\\sqrt{V_i} / \sum\_{j \in i}
  W_j(\hat\gamma_i)\\.

These are the statistics of
[`test.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/test.logistic_fe.md)
on the person-period data. The exact and likelihood-ratio statistics are
not offered: a subject leaves at its event, so a provider's event count
is not a sum of independent trials. Under the model the person-periods
of a subject are uncorrelated, so there is no stabilized variance.

The limits of the score statistic invert the test as in
[`test.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/test.logistic_fe.md);
a provider with no events has no finite lower limit.

## Examples

``` r
s <- sim_kecc(n_providers = 20, provider_size = c(20, 40), rows_per_patient = c(1, 1),
              type = "survival", seed = 1)
s$month <- pmin(ceiling(s$stop * 12), 12)
s$event[s$stop * 12 > 12] <- 0
pd <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider", patient = "patient",
              time = "month", event = "event", message = FALSE)
fit <- fit_discrete_fe(pd)
test(fit)
#> Provider tests: score statistic, naive variance (modified); null median = -3.013; two.sided, level 0.95
#>   20 providers: 3 above, 1 below, 16 not different, 0 not tested
#>  provider estimate se null_value statistic       p flag  lower  upper observed
#>      P001   -3.082 NA     -3.013   -0.2111 0.83277    0 -3.712 -2.452       10
#>      P002   -3.885 NA     -3.013   -1.7819 0.07477    0 -4.840 -2.931        4
#>      P003   -2.535 NA     -3.013    1.5398 0.12362    0 -3.140 -1.931       11
#>      P004   -2.601 NA     -3.013    1.1993 0.23043    0 -3.266 -1.936        9
#>      P005   -3.713 NA     -3.013   -1.7256 0.08442    0 -4.503 -2.922        6
#>      P006   -3.640 NA     -3.013   -1.6634 0.09624    0 -4.376 -2.905        7
#>  expected variance
#>    10.666    9.948
#>     9.299    8.845
#>     7.016    6.695
#>     6.108    5.814
#>    11.758   11.134
#>    12.791   12.123
#>   ... 14 more rows
test(fit, statistic = "wald", variance = "robust")
#> Provider tests: wald statistic, robust variance; null median = -3.013; two.sided, level 0.95
#>   20 providers: 3 above, 1 below, 16 not different, 0 not tested
#>  provider estimate     se null_value statistic       p flag  lower  upper
#>      P001   -3.082 0.3040     -3.013   -0.2269 0.82050    0 -3.677 -2.486
#>      P002   -3.885 0.4948     -3.013   -1.7639 0.07774    0 -4.855 -2.916
#>      P003   -2.535 0.2830     -3.013    1.6865 0.09171    0 -3.090 -1.981
#>      P004   -2.601 0.3123     -3.013    1.3181 0.18748    0 -3.213 -1.989
#>      P005   -3.713 0.4144     -3.013   -1.6895 0.09112    0 -4.525 -2.901
#>      P006   -3.640 0.3824     -3.013   -1.6418 0.10064    0 -4.390 -2.891
#>   ... 14 more rows
```
