# Risk-adjustment models

The risk-adjustment model is where the provider effects come from. This
page explains the choices the package offers (fixed, random or
correlated random provider effects; logistic, linear or Cox models; the
Firth penalty) and when each is the right one. Every fitter takes a
[`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
object and returns a fit that
[`summary()`](https://rdrr.io/r/base/summary.html),
[`standardize()`](https://um-kevinhe.github.io/pprof/reference/standardize.md),
[`test()`](https://um-kevinhe.github.io/pprof/reference/test.md),
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) understand.

## Provider effects as fixed effects

The model is \\\eta\_{ijk} = \gamma_i + Z\_{ijk}^\top\beta\\ with one
free parameter \\\gamma_i\\ per provider. Two things follow.

**The effects are never shrunk.** A random-effect model pulls every
provider toward the average, by more for small providers; a fixed-effect
model reports what the data say about each provider and lets the test
express how uncertain that is. For a measure that reports the provider,
the fixed-effect estimate is the estimand.

**The pooled model is biased.** Dropping the provider effects and
fitting \\Z^\top\beta\\ alone estimates not \\\beta\\ but, in the linear
case, \\\beta + \rho\\\sigma\_\gamma / \sigma_Z\\ when a risk factor is
correlated at \\\rho\\ with the provider effect: a bias that does not
vanish with the sample size. Case mix is related to the provider
(geography, referral patterns, the patients’ circumstances), so \\\rho
\ne 0\\ is the normal state of affairs. For binary outcomes the bias
need not vanish even at \\\rho = 0\\, because of the non-collapsibility
of the logit link.

With thousands of providers the fixed-effect model has thousands of
parameters, which general-purpose software handles by expanding every
provider into an indicator column. The package does not: the block of
the information matrix that belongs to the provider effects is diagonal,
so each effect is updated from its own provider’s records (the serial
block-wise inversion Newton algorithm, SerBIN, of Wu et al. 2022), and a
fit with ten thousand providers takes seconds.

## Binary and binomial outcomes: `fit_logistic_fe()`

``` r

library(pprof)
d <- sim_kecc(n_providers = 60, provider_size = c(15, 60), seed = 7)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_fe(pd)
fit
#> Logistic model with fixed provider effects (SerBIN)
#>   7806 rows, 60 providers, 3 covariates; binary outcome
#>   log-likelihood -4018.33; iterations 5
#>   provider effects: median -0.9008, range -2.48803 to  0.01372
str(coef(fit))
#> List of 2
#>  $ gamma: Named num [1:60] -2.488 -0.341 -0.814 -0.546 -1.687 ...
#>   ..- attr(*, "names")= chr [1:60] "P001" "P002" "P003" "P004" ...
#>  $ beta : Named num [1:3] 0.513 -0.434 0.278
#>   ..- attr(*, "names")= chr [1:3] "x1" "x2" "x3"
```

- `algorithm = "SerBIN"` (default) or `"BAN"` (block ascent Newton, the
  older alternating algorithm; slower, kept for reference).
- `bound`: provider effects are kept within `median(gamma) ± bound` on
  the logit scale (default 10), which is where providers with no events
  or with all events end up.
- A binomial outcome (events out of `n` trials per row, for instance
  hospital days out of days at risk in a month) uses the same fitter
  once
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  was given `n`; every later step accounts for the trials.

The risk-factor coefficients with their tests:

``` r

summary(fit)
#>   covariate   estimate         se statistic            p      lower      upper
#> 1        x1  0.5125504 0.02976985 17.217101 0.000000e+00  0.4542026  0.5708983
#> 2        x2 -0.4343510 0.05811941 -7.473425 7.813411e-14 -0.5482630 -0.3204391
#> 3        x3  0.2779607 0.02810073  9.891581 0.000000e+00  0.2228842  0.3330371
summary(fit, test = "robust_wald")
#>   covariate   estimate         se statistic            p      lower      upper
#> 1        x1  0.5125504 0.03064038 16.727939 0.000000e+00  0.4524964  0.5726045
#> 2        x2 -0.4343510 0.06046739 -7.183227 6.808491e-13 -0.5528649 -0.3158371
#> 3        x3  0.2779607 0.02872757  9.675745 0.000000e+00  0.2216557  0.3342657
```

[`summary()`](https://rdrr.io/r/base/summary.html) offers the
model-based Wald test, the patient-clustered robust Wald test, the
likelihood-ratio test (a refit without the covariate) and the score
test. The variance of the provider effects is `vcov(fit)`’s attribute
`var_gamma`, and
[`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md)
returns the variance the Wald test uses (`"case_mix"`, the variance of
the provider’s fitted average linear predictor, or `"full"`, the
variance of \\\hat\gamma_i\\ alone; see [Tests of provider
effects](https://um-kevinhe.github.io/pprof/articles/testing.md)).

### Providers whose outcomes are all the same: the Firth penalty

A provider with no events has no finite maximum-likelihood effect (it
sits at the lower bound). The Firth penalty, half the log determinant of
the information, gives every provider a finite estimate and reduces the
small-sample bias of the others (Firth 1993). The fit is a `logistic_fe`
object like any other, so every test applies; the likelihood-ratio tests
compare penalized likelihoods.

``` r

firth <- fit_logistic_fe(pd, penalty = "firth")
range(coef(fit)$gamma)
#> [1] -2.48802597  0.01372286
range(coef(firth)$gamma)
#> [1] -2.45472515  0.01455386
```

The kernel follows ParaFLR (Shao and He): provider-level work on
independent threads (`threads`), with the reductions arranged so that
the result does not depend on the thread count.

## Continuous outcomes: `fit_linear_fe()`

A continuous outcome (a laboratory value, a length of stay) uses the
linear fixed-effect model, fitted as a within-provider regression (every
variable centred by provider), which never forms the provider indicators
either.

``` r

dc <- sim_kecc(n_providers = 40, type = "continuous", seed = 8)
pc <- pp_data(dc, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
lin <- fit_linear_fe(pc)
lin
#> Linear model with fixed provider effects
#>   4808 rows, 40 providers, 3 covariates
#>   residual standard deviation 1.114 on 4765 degrees of freedom
#>   provider effects: median 0.06424, range -1.663 to  1.233
summary(lin)
#>    covariate   estimate         se statistic             p      lower
#> x1        x1  0.5484809 0.01651588  33.20930 1.018112e-217  0.5161022
#> x2        x2 -0.3903600 0.03307316 -11.80292  0.000000e+00 -0.4551986
#> x3        x3  0.2955670 0.01614206  18.31036  1.981361e-72  0.2639211
#>         upper
#> x1  0.5808597
#> x2 -0.3255213
#> x3  0.3272129
```

The model-based tests of a linear fit use the \\t\\ distribution with
\\n - m - p\\ degrees of freedom; the standardized measures are
differences rather than ratios.

## Random and correlated random effects

For comparison with the fixed-effect fit, and for the settings where a
shrunken estimate is wanted, the package wraps lme4:

- [`fit_logistic_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md),
  [`fit_linear_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md):
  a random provider intercept, \\\gamma_i \sim N(0, \sigma^2)\\; the
  reported effects are the predicted (shrunken) values,
  `coefficients$alpha`, with conditional variances `post_var`.
- [`fit_logistic_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md),
  [`fit_linear_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md):
  correlated random effects (Mundlak), which add every risk factor’s
  provider mean to the fixed part (`<x>_between`) next to its
  within-provider deviation (`<x>_within`), so that the random effect is
  no longer assumed independent of the case mix. The within coefficients
  equal the fixed-effect coefficients in the linear case.

``` r

re <- fit_logistic_re(pd)
re
#> Logistic model with random provider effects
#>   7806 rows, 60 providers, 3 fixed-effect terms
#>   provider variance 0.2893
#>   log-likelihood -4107.32
cre <- fit_logistic_cre(pd)
summary(cre)
#>               covariate   estimate         se  statistic            p
#> (Intercept) (Intercept) -1.2762688 0.32942663 -3.8742127 1.069700e-04
#> x1_within     x1_within  0.5076747 0.02962586 17.1361988 7.969324e-66
#> x2_within     x2_within -0.4306189 0.05790826 -7.4362256 1.036948e-13
#> x3_within     x3_within  0.2752906 0.02798701  9.8363682 7.849335e-23
#> x1_between   x1_between  0.4545776 0.44587236  1.0195240 3.079543e-01
#> x2_between   x2_between  0.1581450 0.81435872  0.1941958 8.460226e-01
#> x3_between   x3_between -1.1958600 0.83841834 -1.4263286 1.537735e-01
#>                  lower      upper
#> (Intercept) -1.9219331 -0.6306045
#> x1_within    0.4496091  0.5657403
#> x2_within   -0.5441170 -0.3171208
#> x3_within    0.2204370  0.3301441
#> x1_between  -0.4193162  1.3284713
#> x2_between  -1.4379688  1.7542588
#> x3_between  -2.8391298  0.4474097
```

[`standardize()`](https://um-kevinhe.github.io/pprof/reference/standardize.md)
of a random-effect fit offers `numerator = "observed"` or `"predicted"`
(the predicted-to-expected ratio), and
[`test()`](https://um-kevinhe.github.io/pprof/reference/test.md) is the
Wald test of the predicted effect.

## Time-to-event outcomes: `fit_cox()`

For a survival outcome (death, a first admission) or recurrent events
counted on patient-records (admissions, transfusions) the model is a
proportional-hazards model, and the provider effect is handled by
**stratification**: each provider has its own baseline hazard, which
absorbs its effect, and the risk-factor coefficients are estimated
within providers. The national baseline is then the Breslow estimator
over every record with the fitted linear predictor as an offset (the
two-stage model of He and Schaubel 2014). The expected number of events
of every record is \\E = e^{\eta} \[\Lambda_0(\text{stop}) -
\Lambda_0(\text{start})\]\\, and the provider’s measure is its
observed-to-expected ratio.

``` r

s <- sim_kecc(n_providers = 60, provider_size = c(15, 60), type = "survival", seed = 9)
ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider", patient = "patient",
              time = "stop", start = "start", event = "event", cutoff = 0, message = FALSE)
cox <- fit_cox(ps)
cox
#> Cox model stratified by provider (two-stage)
#>   8081 rows, 60 providers, 3 covariates; 1102 events; ties breslow
#>   expected events sum to 1102
summary(cox, variance = "robust")
#>   covariate   estimate         se statistic            p      lower      upper
#> 1        x1  0.4811966 0.03341082 14.402419 4.996364e-47  0.4157126  0.5466806
#> 2        x2 -0.3061843 0.06359755 -4.814404 1.476399e-06 -0.4308332 -0.1815354
#> 3        x3  0.2871641 0.03041721  9.440843 3.697910e-21  0.2275475  0.3467808
head(standardize(cox))
#>   provider observed  expected      ratio person_time
#> 1     P001        7 11.198024 0.62511029    17.76226
#> 2     P002        8 22.315135 0.35850108    41.75757
#> 3     P003        9 14.113514 0.63768669    29.99035
#> 4     P004       11  6.652873 1.65342101    10.98683
#> 5     P005       10 13.595800 0.73552128    24.77177
#> 6     P006        2 31.797469 0.06289809    54.36686
```

- `strata = "provider"` (default) or `"none"` for the pooled baseline.
- `x2`: covariates fitted in a second stage with the first-stage linear
  predictor as offset (calendar-year indicators, say).
- `ties = "breslow"` (default) or `"efron"`.
- Case `weights` and an `offset` declared in
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  are used (recurrent-event measures weight each record by its number of
  events and carry a log-exposure offset).
- The patient column is the cluster of the robust covariance.
- `cutoff = 0`: every provider enters the stratified fit; the providers
  to report are chosen by `test(cox, min_expected = 3)`.

The discrete-time hazard model (one logistic hazard per period, with
provider effects) is fitted by `var_select(ps, model = "discrete")`, see
[Variable
selection](https://um-kevinhe.github.io/pprof/articles/selection.md).

## Which model

| Outcome | Fitter | The provider effect is | Measure |
|----|----|----|----|
| 0/1 per row, or events out of trials | [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md) | a log-odds | ratio or rate |
| continuous | [`fit_linear_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_linear_fe.md) | a difference of means | difference |
| time to event, recurrent events | [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md) | absorbed in the provider’s baseline | observed / expected |
| any, when a shrunken estimate is wanted | `fit_*_re()`, `fit_*_cre()` | a predicted random effect | ratio, rate or difference |

## References

He K, Kalbfleisch JD, Li Y, Li Y (2013). Evaluating hospital readmission
rates in dialysis facilities; adjusting for hospital effects. *Lifetime
Data Analysis*, 19, 490-512.

Wu W, Yang Y, Kang J, He K (2022). Improving large-scale estimation and
inference for profiling health care providers. *Statistics in Medicine*,
41(15), 2840-2853.

Firth D (1993). Bias reduction of maximum likelihood estimates.
*Biometrika*, 80(1), 27-38.

He K, Schaubel DE (2014). Standardized mortality ratio for evaluating
center-specific mortality: assessment and alternative. *Statistics in
Biosciences*, 6(2), 1-19.
