# Standardized measures

A standardized measure summarises the fitted model as one evaluation of
each provider: how its outcomes compare with those expected under a
national norm.
[`standardize()`](https://um-kevinhe.github.io/pprof/reference/standardize.md)
computes the measures of every model class in the same form;
[`test_measure()`](https://um-kevinhe.github.io/pprof/reference/test_measure.md)
tests a measure against a reference value on its own scale.

## Indirect and direct standardization

Let \\O_i\\ be provider \\i\\’s observed total (its outcomes summed, or
its number of events) and \\E_i(a)\\ the total its own patients would
have had at a provider with effect \\a\\, given the fitted risk-factor
coefficients: a sum of fitted means for a logistic or linear model, the
fitted intensity accumulated over the periods at risk for a Cox model.
The norm \\\hat\mu\\ is the median of the provider effects by default
(robust to the outlying providers that are the object of the exercise),
or their size-weighted mean, or a number.

**Indirect** standardization compares the provider’s observed total with
the total its own patients would have had at the norm:

\\\mathrm{ISR}\_i = \frac{O_i}{E_i(\hat\mu)}, \qquad \mathrm{ISD}\_i =
\frac{O_i - E_i(\hat\mu)}{m_i} .\\

The question is “did this provider’s patients do better or worse than
they would have at a typical provider?”. The quantity tested is the
provider effect itself.

**Direct** standardization applies the provider’s fitted effect to the
patients of every provider, and compares the result with the norm
applied to the same patients:

\\\mathrm{DSR}\_i = \frac{\sum\_{i'} E\_{i'}(\hat\gamma_i)}{\sum\_{i'}
E\_{i'}(\hat\mu)}, \qquad \mathrm{DSD}\_i = \frac{\sum\_{i'}
E\_{i'}(\hat\gamma_i) - \sum\_{i'} E\_{i'}(\hat\mu)}{n} .\\

The numerator divided by the total number of observations is the
**directly standardized rate** of provider \\i\\: the rate the whole
population would have had at this provider. Two providers’ direct rates
are comparable without further case-mix adjustment, which is why several
measures report them. The estimand is then a population-average rate, a
smooth monotone function of the provider effect, and testing it needs
the transport described below.

Ratios are used for counts and times to event, differences for
continuous outcomes. A ratio above 1 or a difference above 0 means more
of the outcome than expected.

## `standardize()`

``` r

library(pprof)
d <- sim_kecc(n_providers = 60, provider_size = c(15, 60), seed = 7)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_fe(pd)
ind <- standardize(fit, type = "indirect")
head(ind)
#>   provider observed expected     ratio       rate variance_expected
#> 1     P001       14 50.59789 0.2766914 0.07078564         34.629399
#> 2     P002       44 30.62439 1.4367633 0.36756549         21.551598
#> 3     P003       47 44.45907 1.0571522 0.27045002         29.016328
#> 4     P004       54 42.98972 1.2561143 0.32135027         29.179177
#> 5     P005       35 61.58077 0.5683593 0.14540270         39.925531
#> 6     P006       16 12.27364 1.3036071 0.33350031          8.507391
dir <- standardize(fit, type = "direct")
head(dir)
#>   provider observed  expected     ratio       rate variance_expected
#> 1     P001       14  597.7141 0.2993060 0.07657111          536.3201
#> 2     P002       44 3025.0487 1.5147965 0.38752866         1712.4896
#> 3     P003       47 2261.5853 1.1324914 0.28972397         1497.1159
#> 4     P004       54 2681.7483 1.3428885 0.34354961         1632.1467
#> 5     P005       35 1186.9101 0.5943466 0.15205100          959.0992
#> 6     P006       16 2766.2424 1.3851990 0.35437387         1654.4001
```

The columns: `observed`, `expected`, `ratio`, `rate` (the ratio times
the population rate, on the 0-1 scale) and `variance_expected`
(indirect: \\\sum N p_0(1 - p_0)\\, the variance of the expected count
under the null; direct: the delta-method numerator). Attributes record
the type, the null value, the population size and rate.

- `null`: `"median"` (default), `"mean"` (weighted by rows) or a number.
- `n_pop`: the denominator of the rates; by default the number of trials
  among the fitted rows.
- `providers`: a subset to report; the norm and the expectations always
  use every fitted row.

### Confidence limits of the measures

Limits of the ratio and the rate are mapped from limits of the provider
effect, obtained by inverting one of the tests (`ci = "exact"`,
`"score"`, `"wald"` or `"robust_wald"`):

``` r

head(standardize(fit, type = "indirect", level = 0.95, ci = "score"))
#>   provider observed expected     ratio       rate variance_expected lower_ratio
#> 1     P001       14 50.59789 0.2766914 0.07078564         34.629399   0.1673478
#> 2     P002       44 30.62439 1.4367633 0.36756549         21.551598   1.1289012
#> 3     P003       47 44.45907 1.0571522 0.27045002         29.016328   0.8321438
#> 4     P004       54 42.98972 1.2561143 0.32135027         29.179177   1.0091884
#> 5     P005       35 61.58077 0.5683593 0.14540270         39.925531   0.4204167
#> 6     P006       16 12.27364 1.3036071 0.33350031          8.507391   0.8598169
#>   upper_ratio lower_rate upper_rate
#> 1   0.4463502 0.04281239  0.1141893
#> 2   1.7717141 0.28880551  0.4532556
#> 3   1.3092334 0.21288640  0.3349397
#> 4   1.5259656 0.25817952  0.3903860
#> 5   0.7526987 0.10755471  0.1925620
#> 6   1.8270113 0.21996598  0.4674022
```

### The standard error of a direct rate

A direct rate is a function of one provider effect, so any standard
error of \\\hat\gamma_i\\ transports to it by the delta method:
`se_rate = numer_se / n_pop * se(gamma_i)`, with \\\hat\beta\\ held
fixed. Give the standard errors of a Wald test as `gamma_se`:

``` r

tst <- test(fit, statistic = "wald", variance = "robust")
dir <- standardize(fit, type = "direct", gamma_se = stats::setNames(tst$se, tst$provider))
head(dir[, c("provider", "rate", "se_rate")])
#>   provider       rate    se_rate
#> 1     P001 0.07657111 0.02144172
#> 2     P002 0.38752866 0.05324641
#> 3     P003 0.28972397 0.04118237
#> 4     P004 0.34354961 0.03785637
#> 5     P005 0.15205100 0.02336948
#> 6     P006 0.35437387 0.08554653
```

## Testing a direct measure: `test_measure()`

The measures that report direct rates test the rate against the national
rate on the rate scale (or on the logit scale), with the transported
standard error.
[`test_measure()`](https://um-kevinhe.github.io/pprof/reference/test_measure.md)
is the measure-level z test \\z_i = \[g(M_i) - g(M\_{ref})\] /
\[g'(M_i)\\ SE_i\]\\ with \\g\\ the identity, the logit or the log:

``` r

mt <- test_measure(estimate = dir$rate, se = dir$se_rate,
                   reference = attr(dir, "population_rate"), transform = "logit",
                   providers = dir$provider)
head(mt)
#>   provider   estimate         se reference estimate_t      se_t reference_t
#> 1     P001 0.07657111 0.02144172 0.2558288 -2.4898739 0.3032432   -1.067762
#> 2     P002 0.38752866 0.05324641 0.2558288 -0.4577124 0.2243369   -1.067762
#> 3     P003 0.28972397 0.04118237 0.2558288 -0.8967250 0.2001243   -1.067762
#> 4     P004 0.34354961 0.03785637 0.2558288 -0.6475156 0.1678601   -1.067762
#> 5     P005 0.15205100 0.02336948 0.2558288 -1.7186045 0.1812550   -1.067762
#> 6     P006 0.35437387 0.08554653 0.2558288 -0.5998681 0.3739036   -1.067762
#>            z            p flag      lower     upper    lower_t     upper_t
#> 1 -4.6896726 2.736426e-06   -1 0.04376289 0.1306154 -3.0842198 -1.89552808
#> 2  2.7193474 6.541087e-03    1 0.28958414 0.4954951 -0.8974046 -0.01802007
#> 3  0.8546557 3.927418e-01    0 0.21602866 0.3764864 -1.2889614 -0.50448868
#> 4  2.5035532 1.229532e-02    1 0.27358375 0.4210375 -0.9765154 -0.31851575
#> 5 -3.5907544 3.297222e-04   -1 0.11166379 0.2036962 -2.0738578 -1.36335129
#> 6  1.2513769 2.107970e-01    0 0.20871218 0.5331935 -1.3327056  0.13296939
table(mt$flag)
#> 
#> -1  0  1 
#> 17 28 15
```

Its null distribution is \\N(\text{location}, \text{scale}^2)\\: the
theoretical null by default, a known inflation such as `scale = 1.8`, or
an empirical null through `null_model` (see [Reference
distributions](https://um-kevinhe.github.io/pprof/articles/null.md)).
The intervals invert the test and are mapped back to the measure scale.

## Continuous outcomes

For a linear fit the measures are differences: indirect \\(O_i -
E_i)/m_i\\, direct relative to the population mean of the effects.

``` r

dc <- sim_kecc(n_providers = 40, type = "continuous", seed = 8)
pc <- pp_data(dc, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
head(standardize(fit_linear_fe(pc)))
#>   provider    observed  expected difference
#> 1     P001  -57.119267 -28.33962 -0.1998587
#> 2     P002    7.166393 -13.64698  0.1118998
#> 3     P003   11.282420 -23.39783  0.3181675
#> 4     P004   56.397734 -24.26554  0.9379450
#> 5     P005 -109.383890 -11.77728 -0.5483518
#> 6     P006   85.828761 -16.77811  0.5965516
```

## Survival outcomes

For a Cox fit the indirect measure is the observed-to-expected ratio
\\O_j / E_j\\ with \\E_j\\ the events expected at the national baseline
for the provider’s own records, together with the person-time (the
grouping variable of the survival measures’ empirical null); the direct
measure is \\E^{(j)} / O\\ with \\E^{(j)}\\ the events expected if every
patient had provider \\j\\’s own baseline, computed from the risk-set
sums.

``` r

s <- sim_kecc(n_providers = 60, provider_size = c(15, 60), type = "survival", seed = 9)
ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider", patient = "patient",
              time = "stop", start = "start", event = "event", cutoff = 0, message = FALSE)
cox <- fit_cox(ps)
head(standardize(cox))
#>   provider observed  expected      ratio person_time
#> 1     P001        7 11.198024 0.62511029    17.76226
#> 2     P002        8 22.315135 0.35850108    41.75757
#> 3     P003        9 14.113514 0.63768669    29.99035
#> 4     P004       11  6.652873 1.65342101    10.98683
#> 5     P005       10 13.595800 0.73552128    24.77177
#> 6     P006        2 31.797469 0.06289809    54.36686
head(standardize(cox, type = "direct"))
#>   provider observed   expected      ratio n_pop
#> 1     P001     1102  666.21886 0.60455432  8081
#> 2     P002     1102  402.91185 0.36561873  8081
#> 3     P003     1102  776.90991 0.70499992  8081
#> 4     P004     1102 1648.87585 1.49625758  8081
#> 5     P005     1102  731.20286 0.66352346  8081
#> 6     P006     1102   62.61549 0.05681986  8081
```

The ratios of a Cox fit are tested by `test(cox)`, which treats the
observed count as Poisson with the expected count as its mean (see
[Tests of provider
effects](https://um-kevinhe.github.io/pprof/articles/testing.md)).

## An arbitrary standard population

[`direct_expectation()`](https://um-kevinhe.github.io/pprof/reference/direct_expectation.md)
evaluates the direct numerator for any provider effects and any standard
population (a linear predictor and trials per row), for instance a
national population from another year:

``` r

direct_expectation(gamma = coef(fit)$gamma[1:3], linear_predictor = fit$linear_predictor)
#>    expected  variance
#> 1  597.7141  536.3201
#> 2 3025.0487 1712.4896
#> 3 2261.5853 1497.1159
```

## Reference

He K (2019). Indirect and direct standardization for evaluating
transplant centers. *Journal of Hospital Administration*, 8(1), 9-14.
