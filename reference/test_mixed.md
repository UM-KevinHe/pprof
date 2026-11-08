# Tests of the provider effects of a random-effect fit

Tests \\H_0: \alpha_i = \alpha_0\\ for every provider (or a subset) with
the statistic \\(\hat\alpha_i - \alpha_0) / \sqrt{v_i}\\, where
\\\hat\alpha_i\\ is the predicted random effect and \\v_i\\ its
conditional variance given the data (`post_var`), referred to the
standard normal.

## Usage

``` r
# S3 method for class 'logistic_re'
test(
  fit,
  null = "median",
  providers = NULL,
  level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  null_model = NULL,
  ...
)

# S3 method for class 'linear_re'
test(
  fit,
  null = "median",
  providers = NULL,
  level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  null_model = NULL,
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

- null:

  the null value \\\alpha_0\\: `"median"` (default) or `"mean"`
  (weighted by rows) of the predicted effects, or a number (0 is the
  population intercept).

- providers:

  optional subset of provider identifiers to report.

- level:

  confidence level; the significance level is `1 - level`.

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

- null_model:

  optional null model (a `pp_null` object or a function of the result
  returning one), applied through
  [`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md).

- ...:

  unused.

## Value

A data frame of class `pp_test` with the columns described in
[`test.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/test.logistic_fe.md);
`estimate` is \\\hat\alpha_i\\.

## Examples

``` r
# \donttest{
d <- sim_kecc(n_providers = 20, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
test(fit_logistic_re(pd))
#> Provider tests: wald statistic, conditional variance; null median = -0.05977; two.sided, level 0.95
#>   20 providers: 5 above, 3 below, 12 not different, 0 not tested
#>  provider estimate     se null_value statistic         p flag   lower    upper
#>      P001 -0.40407 0.2728   -0.05977   -1.2621 2.069e-01    0 -0.9388  0.13062
#>      P002 -0.44762 0.1847   -0.05977   -2.0995 3.577e-02   -1 -0.8097 -0.08555
#>      P003  0.72599 0.2709   -0.05977    2.9000 3.731e-03    1  0.1949  1.25703
#>      P004  0.63557 0.1616   -0.05977    4.3019 1.694e-05    1  0.3188  0.95237
#>      P005  0.06144 0.1940   -0.05977    0.6247 5.322e-01    0 -0.3189  0.44175
#>      P006 -0.31541 0.1675   -0.05977   -1.5264 1.269e-01    0 -0.6437  0.01285
#>   ... 14 more rows
# }
```
