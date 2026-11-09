# Calibrate a test result with a null model

Refers the statistics of a
[`test()`](https://um-kevinhe.github.io/pprof/reference/test.md) or
[`test_measure()`](https://um-kevinhe.github.io/pprof/reference/test_measure.md)
result to a null model: the adjusted statistic \\(z -
\mathrm{null\\mean}) / \mathrm{null\\sd}\\ replaces the standard normal
in the p-values and flags, and the confidence limits are those of the
inverted calibrated test, \\\hat\theta - (\mathrm{null\\mean} \pm c\\
\mathrm{null\\sd})\\ se\\, so that a limit excludes the null value
exactly when the provider is flagged.

## Usage

``` r
calibrate(
  x,
  null_model,
  level = NULL,
  alternative = NULL,
  critical = NULL,
  ci_location = NULL
)
```

## Arguments

- x:

  a test result.

- null_model:

  a `pp_null` object, or a function of `x` returning one (for instance
  `function(t) null_model("empirical", t, n_groups = 4)`).

- level, alternative:

  default to the settings of `x`.

- critical:

  optional critical value for the flags and limits instead of the normal
  quantile for `level`.

- ci_location:

  location at which the confidence limits are centred: `NULL` (default)
  uses the null mean, so that a limit excludes the null value exactly
  when the provider is flagged; a number or one value per provider
  centres the limits elsewhere, for instance `0` for limits of the form
  \\\hat\theta \pm c\\ \mathrm{null\\sd}\\ se\\ that ignore the null
  location.

## Value

`x` with the columns `z_adjusted`, `null_mean`, `null_sd` and
`null_group` added, `p` and `flag` recomputed from the adjusted
statistic, and `lower`/`upper` recomputed where a standard error is
available (for
[`test_measure()`](https://um-kevinhe.github.io/pprof/reference/test_measure.md)
results on the transformed scale, mapped back; for the log-ratio Wald
and mid-p tests of
[`test.cox()`](https://um-kevinhe.github.io/pprof/reference/test.cox.md)
on the ratio scale, the mid-p limits by root finding); the settings are
recorded in the `test` attribute. The exact Poisson test cannot be
calibrated.

## Details

For a Wald-type test the limits are \\\hat\theta - (\mathrm{null\\mean}
\pm c\\ \mathrm{null\\sd})\\ se\\; for a modified score test of
[`test.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/test.logistic_fe.md)
they invert the calibrated statistic (the set of null values at which
the adjusted statistic lies within \\\pm c\\), for a mid-p test the
calibrated mid-p rule; the exact Poisson test cannot be calibrated.

## Examples

``` r
d <- sim_kecc(n_providers = 40, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
tst <- test(fit_logistic_fe(pd), statistic = "wald", variance = "robust")
calibrate(tst, null_model("fixed", sd = 1.5))
#> Provider tests: wald statistic, robust variance; null median = -0.965; two.sided, level 0.95
#>   40 providers: 6 above, 4 below, 30 not different, 0 not tested
#>  provider estimate     se null_value statistic z_adjusted null_mean null_sd
#>      P001  -0.7418 0.4278     -0.965    0.5215     0.3477         0     1.5
#>      P002  -0.3816 0.1720     -0.965    3.3920     2.2614         0     1.5
#>      P003  -1.7131 0.3869     -0.965   -1.9338    -1.2892         0     1.5
#>      P004  -1.9726 0.2968     -0.965   -3.3953    -2.2636         0     1.5
#>      P005  -1.5766 0.2120     -0.965   -2.8847    -1.9231         0     1.5
#>      P006  -0.6088 0.1739     -0.965    2.0480     1.3653         0     1.5
#>  null_group       p flag   lower    upper
#>          NA 0.72807    0 -1.9996  0.51589
#>          NA 0.02374    1 -0.8872  0.12396
#>          NA 0.19732    0 -2.8504 -0.57575
#>          NA 0.02360   -1 -2.8452 -1.10011
#>          NA 0.05446    0 -2.2000 -0.95325
#>          NA 0.17214    0 -1.1200 -0.09763
#>   ... 34 more rows
```
