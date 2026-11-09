# Tests of provider effects

Profiling reduces to one hypothesis test per provider, of the null that
the provider’s effect equals a reference effect \\\gamma_0\\, the median
of the provider effects. Three difficulties shape the test: small
providers, some with outcomes that are all identical; observations
correlated within patient; and a reference distribution that may not
describe how providers ordinarily vary. They correspond to three
choices, made independently: the **statistic**, its **variance** and the
**reference distribution**. This page covers the first two; the third is
[Reference
distributions](https://um-kevinhe.github.io/pprof/articles/null.md).

## The taxonomy

|  | naive (working independence) | robust (patient-clustered sandwich) | stabilized robust (pooled correlation) |
|----|----|----|----|
| Wald | \\Z_W\\ | \\Z\_{RW}\\ | \\Z\_{SRW}\\ |
| score | \\Z_S\\ | \\Z\_{RS}\\ | \\Z\_{SRS}\\ |
| exact | \\Z_E\\ |  |  |

**Three statistics.** A *Wald* statistic compares the estimated effect
with the reference, scaled by its standard error; it needs the estimate
to exist, and at a provider whose outcomes are all identical it does
not. A *score* statistic evaluates the estimating function at the null,
so nothing is estimated for the provider under test and it is finite
exactly where the Wald statistic breaks. An *exact* test uses the
finite-sample distribution of the provider’s total under the null:
Poisson-binomial for a binary outcome, from which a mid-p value and an
interval follow; it assumes the observations within the provider are
independent.

**Three variances.** The *naive* variance treats a patient’s repeated
observations as independent, which they are not, so it is too small. The
*robust* variance is the patient-clustered sandwich, consistent whatever
the correlation, but at small providers it is variable and biased
downward, so the robust test flags too many small providers. The
*stabilized* variance (Pan 2001) replaces the provider’s own estimate of
the within-patient correlation with one pooled over every patient in the
data with the same number of observations: the correlation is a feature
of the outcome rather than of the provider, so this is a natural thing
to borrow, and the estimator is far less variable at small providers.

## `test()`

``` r

library(pprof)
d <- sim_kecc(n_providers = 80, provider_size = c(12, 60), rows_per_patient = c(1, 8),
              seed = 12)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_fe(pd)
tst <- test(fit, statistic = "score", variance = "stabilized")
tst
#> Provider tests: score statistic, stabilized variance (modified); null median = -0.9478; two.sided, level 0.95
#>   80 providers: 15 above, 20 below, 45 not different, 0 not tested
#>  provider estimate se null_value statistic         p flag lower upper observed
#>      P001  -0.3947 NA    -0.9478    1.7121 8.687e-02    0    NA    NA       21
#>      P002  -0.1281 NA    -0.9478    4.3560 1.325e-05    1    NA    NA       71
#>      P003  -0.3444 NA    -0.9478    2.6875 7.198e-03    1    NA    NA       46
#>      P004  -0.6544 NA    -0.9478    1.4487 1.474e-01    0    NA    NA       57
#>      P005  -1.0141 NA    -0.9478   -0.3859 6.996e-01    0    NA    NA       72
#>      P006  -0.6903 NA    -0.9478    0.8584 3.907e-01    0    NA    NA       24
#>  expected variance
#>     14.77    13.22
#>     42.78    41.96
#>     31.15    30.54
#>     47.31    44.75
#>     75.24    70.42
#>     20.20    19.60
#>   ... 74 more rows
```

Every result has the same columns: `provider`, `estimate` (the provider
effect), `se` (Wald statistics only), `null_value`, `statistic`, `p`,
`flag` (`1` above the null, `-1` below, `0` not different, `NA` not
testable), `lower`, `upper`. Attributes keep the settings and the
provider sizes.

``` r

combos <- expand.grid(statistic = c("exact", "score", "wald"),
                      variance = c("naive", "robust", "stabilized"), stringsAsFactors = FALSE)
combos <- combos[!(combos$statistic == "exact" & combos$variance != "naive"), ]
combos$flagged <- mapply(function(s, v) {
  sum(test(fit, statistic = s, variance = v)$flag != 0, na.rm = TRUE)
}, combos$statistic, combos$variance)
combos
#>   statistic   variance flagged
#> 1     exact      naive      42
#> 2     score      naive      43
#> 3      wald      naive      43
#> 5     score     robust      36
#> 6      wald     robust      42
#> 8     score stabilized      35
#> 9      wald stabilized      40
```

Arguments:

- `null`: `"median"` (default), `"mean"` (weighted by rows) or a number.
- `alternative`: `"two.sided"`, `"greater"` (the provider is above the
  reference) or `"less"`. A measure whose consequential tail is “worse”
  tests one-sided, at the one-sided level.
- `level`: the confidence level; the significance level is `1 - level`.
- `providers`: a subset to report (every provider still enters the
  stabilized pooling).
- `exact_method`: `"poibin"` (the Poisson-binomial distribution) or
  `"bootstrap"`; `n_resample`.
- `score_type`: `"modified"` (the other parameters held at their
  full-model estimates) or `"restricted"` (re-solved under the null).
- `wald_variance`: `"case_mix"` (default) or `"full"`, see below.
- `statistic = "lr"`: the likelihood-ratio statistic of one provider
  effect with the other parameters fixed (ParaFLR), the penalized
  likelihood for a Firth fit.

### The Wald variance

The naive Wald test uses by default the variance of the provider’s
fitted average linear predictor, \\\mathrm{Var}(\hat\gamma_i + \bar
Z_i^\top \hat\beta)\\, which does not depend on where the covariates are
centred; the `"full"` variance of \\\hat\gamma_i\\ alone does, and with
uncentred covariates it changes the flags.

### Confidence limits

[`confint()`](https://rdrr.io/r/stats/confint.html) inverts the chosen
test for the provider effects (`type = "exact"`, `"score"`, `"wald"`,
`"robust_wald"`); `standardize(level = )` maps those limits to the
measures.

``` r

head(confint(fit, type = "score"))
#>   provider   estimate      lower      upper
#> 1     P001 -0.3946980 -0.9508951  0.1617724
#> 2     P002 -0.1280870 -0.4455044  0.1893382
#> 3     P003 -0.3444118 -0.7207547  0.0319795
#> 4     P004 -0.6543873 -0.9869470 -0.3217850
#> 5     P005 -1.0141023 -1.2955286 -0.7326444
#> 6     P006 -0.6903412 -1.1846614 -0.1958079
```

## Which test

The rules the measures follow, derived from the data:

- a provider with all-identical outcomes has no Wald statistic: use the
  score or the exact statistic (or a Firth fit);
- with repeated observations per patient the naive variance is wrong:
  use the robust or the stabilized variance;
- with few patients per provider the robust variance is unstable: use
  the stabilized one;
- when the statistics spread more than \\N(0,1)\\ and the spread grows
  with provider size, calibrate with an empirical null.

The default `test(fit)` is the exact test; the production choice for a
patient-month measure is the robust or stabilized score test (or the
robust Wald test on a direct rate) under an empirical null.

## Continuous and random-effect models

A linear fit has the score and Wald statistics with the naive (a \\t\\
reference with \\n - m - p\\ degrees of freedom) or the robust variance;
a random-effect fit has the Wald test of the predicted effect.

``` r

dc <- sim_kecc(n_providers = 40, type = "continuous", seed = 8)
pc <- pp_data(dc, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
test(fit_linear_fe(pc), variance = "robust")
#> Provider tests: wald statistic, robust variance; null median = 0.06424; two.sided, level 0.95
#>   40 providers: 15 above, 13 below, 12 not different, 0 not tested
#>  provider estimate      se null_value statistic         p flag    lower   upper
#>      P001  -0.1356 0.13634    0.06424    -1.466 1.427e-01    0 -0.40284  0.1316
#>      P002   0.1761 0.10747    0.06424     1.041 2.978e-01    0 -0.03451  0.3868
#>      P003   0.3824 0.10852    0.06424     2.932 3.370e-03    1  0.16970  0.5951
#>      P004   1.0022 0.14592    0.06424     6.428 1.295e-10    1  0.71618  1.2882
#>      P005  -0.4841 0.09295    0.06424    -5.899 3.654e-09   -1 -0.66630 -0.3019
#>      P006   0.6608 0.11217    0.06424     5.318 1.049e-07    1  0.44093  0.8806
#>   ... 34 more rows
```

## Survival outcomes: tests of the observed-to-expected ratio

For a Cox fit the provider’s measure is \\O_j / E_j\\, and
[`test()`](https://um-kevinhe.github.io/pprof/reference/test.md) treats
\\O_j\\ as Poisson with mean \\E_j\\:

- `statistic = "midp"` (default): the two-sided mid-p Poisson test
  turned into a z statistic, whose limits invert the (calibrated) test
  by root finding, so that a limit excludes 1 exactly when the provider
  is flagged;
- `"exact"`: the two-sided Poisson test with the exact chi-square limits
  (\\E_j \< 100\\) or Byar’s approximation; theoretical null only;
- `"wald"`: \\\log(O_j / E_j) / \sqrt{\phi / O_j}\\ with
  `dispersion = "none"` (\\\phi = 1\\) or `"deviance"` (the Poisson
  deviance of the patient-level counts over its degrees of freedom).

`min_expected` drops providers with fewer expected events than the
measure reports (3, typically) before an empirical null is fitted.

``` r

s <- sim_kecc(n_providers = 80, provider_size = c(12, 80), type = "survival", seed = 9)
ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider", patient = "patient",
              time = "stop", start = "start", event = "event", cutoff = 0, message = FALSE)
cox <- fit_cox(ps)
test(cox, min_expected = 3)
#> Provider tests: midp statistic, poisson variance; null ratio = 1; two.sided, level 0.95
#>   80 providers: 14 above, 14 below, 52 not different, 0 not tested
#>  provider estimate se null_value statistic         p flag   lower upper
#>      P001   1.1315 NA          1    0.7526 0.4517142    0 0.80450 1.550
#>      P002   0.4261 NA          1   -3.6620 0.0002503   -1 0.24255 0.698
#>      P003   0.4421 NA          1   -1.7871 0.0739269    0 0.14048 1.066
#>      P004   1.2068 NA          1    1.1480 0.2509728    0 0.86618 1.639
#>      P005   1.1863 NA          1    0.7769 0.4372359    0 0.74497 1.800
#>      P006   0.3083 NA          1   -1.9199 0.0548665    0 0.05169 1.019
#>  observed expected person_time
#>        36   31.815       61.67
#>        14   32.854       65.23
#>         4    9.047       18.54
#>        38   31.489       54.34
#>        20   16.859       30.39
#>         2    6.487       13.71
#>   ... 74 more rows
test(cox, statistic = "wald", dispersion = "deviance", min_expected = 3)
#> Provider tests: wald statistic, poisson variance; null ratio = 1; two.sided, level 0.95
#>   80 providers: 12 above, 10 below, 58 not different, 0 not tested
#>  provider estimate     se null_value statistic        p flag   lower  upper
#>      P001   1.1315 0.1925          1    0.6420 0.520905    0 0.77592 1.6501
#>      P002   0.4261 0.3087          1   -2.7635 0.005718   -1 0.23270 0.7803
#>      P003   0.4421 0.5775          1   -1.4133 0.157560    0 0.14256 1.3712
#>      P004   1.2068 0.1874          1    1.0031 0.315811    0 0.83588 1.7422
#>      P005   1.1863 0.2583          1    0.6615 0.508297    0 0.71510 1.9680
#>      P006   0.3083 0.8167          1   -1.4408 0.149654    0 0.06221 1.5281
#>  observed expected person_time
#>        36   31.815       61.67
#>        14   32.854       65.23
#>         4    9.047       18.54
#>        38   31.489       54.34
#>        20   16.859       30.39
#>         2    6.487       13.71
#>   ... 74 more rows
```

## Tests of the risk-factor coefficients

[`summary()`](https://rdrr.io/r/base/summary.html) tests the
coefficients: Wald with the model-based or the patient-clustered robust
variance, the likelihood ratio and the score test (refits without the
covariate) for the logistic fit; Wald tests of either stage for the Cox
fit; `plot(summary(fit))` is a forest plot.

``` r

summary(fit, test = "robust_wald")
#>   covariate   estimate         se statistic            p      lower      upper
#> 1        x1  0.4951529 0.02386172 20.750933 0.000000e+00  0.4483848  0.5419210
#> 2        x2 -0.4424828 0.04840342 -9.141562 6.155656e-20 -0.5373518 -0.3476139
#> 3        x3  0.2922481 0.02199266 13.288439 0.000000e+00  0.2491433  0.3353530
```

## References

Wu W, Yang Y, Kang J, He K (2022). Improving large-scale estimation and
inference for profiling health care providers. *Statistics in Medicine*,
41(15), 2840-2853.

Pan W (2001). On the robust variance estimator in generalised estimating
equations. *Biometrika*, 88(3), 901-906.

He K, Kalbfleisch JD, Li Y, Li Y (2013). Evaluating hospital readmission
rates in dialysis facilities; adjusting for hospital effects. *Lifetime
Data Analysis*, 19, 490-512.
