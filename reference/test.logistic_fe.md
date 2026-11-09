# Tests of the provider effects of a logistic fixed-effect fit

Tests \\H_0: \gamma_i = \gamma_0\\ for every provider (or a subset) with
a chosen statistic and variance, and flags the providers above or below
the null value.

## Usage

``` r
# S3 method for class 'logistic_fe'
test(
  fit,
  statistic = c("exact", "score", "wald", "lr"),
  variance = c("naive", "robust", "stabilized"),
  null = "median",
  providers = NULL,
  level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  exact_method = c("poibin", "bootstrap"),
  score_type = c("modified", "restricted"),
  wald_variance = c("case_mix", "full"),
  n_resample = 10000,
  null_model = NULL,
  interval = TRUE,
  ...
)
```

## Arguments

- fit:

  a model fitted by
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md).

- statistic:

  `"exact"` (default), `"score"`, `"wald"` or `"lr"` (the
  likelihood-ratio statistic of \\\gamma_i\\ with every other parameter
  held at its full-model estimate; the penalized likelihood of a Firth
  fit, as in ParaFLR).

- variance:

  `"naive"` (default; model-based), `"robust"` (patient-clustered
  sandwich) or `"stabilized"` (sandwich with one within-patient
  correlation matrix pooled over all patients). The clustered variances
  need a patient column in
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  and are not available for the exact test.

- null:

  the null value \\\gamma_0\\: `"median"` (default) or `"mean"`
  (weighted by rows) of the fitted effects, or a number.

- providers:

  optional subset of provider identifiers to report. Every provider
  still enters the stabilized pooling.

- level:

  confidence level; the significance level is `1 - level`.

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

- exact_method:

  `"poibin"` (default; the Poisson-binomial distribution of the
  provider's event count, each row contributing `n` Bernoulli trials) or
  `"bootstrap"` (the same by simulation).

- score_type:

  `"modified"` (default; \\\hat\beta\\ held at the full-model estimate)
  or `"restricted"` (the other parameters re-solved under the null;
  naive variance only).

- wald_variance:

  for `statistic = "wald"` with the naive variance: `"case_mix"`
  (default) or `"full"`, see
  [`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md).

- n_resample:

  number of draws for `exact_method = "bootstrap"`.

- null_model:

  optional null model: a `pp_null` object or a function of the result
  returning one, for instance
  `function(t) null_model("empirical", t, n_groups = 4)`; the result is
  then passed through
  [`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md).

- interval:

  for the score statistics, compute the confidence limits of
  \\\gamma_i\\ by inverting the test (default `TRUE`; `FALSE` saves the
  root finding on very large data). The Wald statistics always carry
  their limits.

- ...:

  unused.

## Value

A data frame of class `pp_test` with one row per tested provider:
`provider`, `estimate` (\\\hat\gamma_i\\), `se` (Wald statistics only),
`null_value`, `statistic` (a z statistic; for the exact test the normal
quantile of the mid-p, for `"lr"` the signed root of the chi-square
statistic), `p`, `flag` (`1` above the null, `-1` below, `0` not
different, `NA` not testable) and, for the Wald and (with
`interval = TRUE`) the modified score statistics, `lower` and `upper`
limits for \\\gamma_i\\. The score and exact statistics also carry
`observed` (the provider's events), `expected` (under the null value)
and `variance` (the variance the statistic is referred to; `NA` for the
restricted score), which the funnel plot draws on. Attributes `test`
(the settings) and `sizes` (rows and patients per provider).

## Details

Row residuals \\r_j = Y_j - N_j p_j\\ and variances \\v_j = N_j p_j(1 -
p_j)\\ are evaluated at the fitted \\\hat p_j\\ for the Wald statistics
and at \\p_j = \mathrm{expit}(\gamma_0 + Z_j^\top\hat\beta)\\ for the
score statistics. With \\U_i = \sum\_{j \in i} r_j\\ and \\I_i =
\sum\_{j \in i} v_j\\:

- naive score: \\U_i / \sqrt{I_i}\\;

- robust score: \\U_i / \sqrt{V_i}\\, \\V_i = \sum\_{k \in i} (\sum\_{j
  \in k} r_j)^2\\ over the provider's patients \\k\\;

- stabilized score: \\V_i\\ replaced by \\\sum\_{k \in i} a_k^\top \hat
  R a_k\\ with \\a_k\\ the \\\sqrt{v_j}\\ of patient \\k\\'s rows and
  \\\hat R\\ pooled over all patients by row position;

- naive Wald: \\(\hat\gamma_i - \gamma_0)/\sqrt{\mathrm{Var}}\\ with the
  variance of
  [`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md);

- robust and stabilized Wald: standard error \\\sqrt{V_i} / I_i\\;

- likelihood ratio: \\-2\[\ell^\*(\hat\gamma\_{(i)}, \hat\beta) -
  \ell^\*(\hat\gamma, \hat\beta)\]\\ with \\\hat\gamma\_{(i)}\\ the
  fitted effects except \\\gamma_i = \gamma_0\\, and \\\ell^\*\\ the
  (penalized, for a Firth fit) log-likelihood; the other parameters are
  not re-estimated under the null. The statistic reported is the signed
  square root, so that two-sided p-values equal the chi-square tail.

The stabilized statistics take a patient's rows in their order in the
data, so sort by time within patient before
[`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md).
A patient with more than 5,000 rows stops the stabilized statistics.

The limits of a score statistic invert the test: with \\z_i(g)\\ the
statistic at the null value \\g\\ (the same residuals and variance,
naive, robust or stabilized, evaluated at \\g\\; the correlation matrix
of the stabilized variance is held at the test's null), the interval is
the set of \\g\\ with \\-c \le z_i(g) \le c\\, so it excludes the null
value exactly when the provider is flagged. A provider with no events
has no finite lower limit, one with all events no finite upper limit,
and a patient-clustered statistic, which is bounded by the square root
of the number of events (non-events), can leave a limit infinite for a
provider with very few of them.
[`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md)
inverts the calibrated statistic the same way.

## References

Wu W, Yang Y, Kang J, He K (2022). Improving large-scale estimation and
inference for profiling health care providers. *Statistics in Medicine*,
41(15), 2840-2853.

## Examples

``` r
d <- sim_kecc(n_providers = 20, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_fe(pd)
test(fit, statistic = "score", variance = "robust")
#> Provider tests: score statistic, robust variance (modified); null median = -1.08; two.sided, level 0.95
#>   20 providers: 5 above, 3 below, 12 not different, 0 not tested
#>  provider  estimate se null_value statistic         p flag   lower   upper
#>      P001 -1.586735 NA      -1.08   -1.7858 0.0741319    0 -2.5245 -1.0209
#>      P002 -1.524872 NA      -1.08   -2.3947 0.0166351   -1 -1.9638 -1.1647
#>      P003  0.006187 NA      -1.08    2.3114 0.0208110    1 -0.8075  0.7553
#>      P004 -0.298864 NA      -1.08    3.7175 0.0002012    1 -0.6257  0.0210
#>      P005 -0.935401 NA      -1.08    0.6637 0.5068676    0 -1.4430 -0.5364
#>      P006 -1.362773 NA      -1.08   -1.5127 0.1303440    0 -1.8155 -0.9975
#>  observed expected variance
#>        11    16.05    7.982
#>        31    43.30   26.392
#>        21    11.12   18.268
#>        63    38.55   43.259
#>        34    30.80   23.271
#>        42    51.42   38.770
#>   ... 14 more rows
```
