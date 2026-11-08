# Tests of the provider effects of a linear fixed-effect fit

Wald tests of \\H_0: \gamma_i = \gamma_0\\ for every provider (or a
subset) with the model-based or the patient-clustered standard error.

## Usage

``` r
# S3 method for class 'linear_fe'
test(
  fit,
  variance = c("naive", "robust"),
  null = "median",
  providers = NULL,
  level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  wald_variance = c("case_mix", "full"),
  null_model = NULL,
  ...
)
```

## Arguments

- fit:

  a model fitted by
  [`fit_linear_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_linear_fe.md).

- variance:

  `"naive"` (default; model-based, see
  [`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md))
  or `"robust"` (patient-clustered sandwich with \\\hat\beta\\ held
  fixed; needs a patient column in
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)).

- null:

  the null value \\\gamma_0\\: `"median"` (default) or `"mean"`
  (weighted by rows) of the fitted effects, or a number.

- providers:

  optional subset of provider identifiers to report.

- level:

  confidence level; the significance level is `1 - level`.

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

- wald_variance:

  for `variance = "naive"`: `"case_mix"` (default) or `"full"`, see
  [`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md).

- null_model:

  optional null model (a `pp_null` object or a function of the result
  returning one), applied through
  [`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md).

- ...:

  unused.

## Value

A data frame of class `pp_test` with the columns described in
[`test.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/test.logistic_fe.md).

## Details

The naive statistic \\(\hat\gamma_i - \gamma_0) / se_i\\ is referred to
the t distribution with \\n - m - p\\ degrees of freedom, the robust
statistic to the standard normal. The robust standard error is
\\\sqrt{\sum\_{k \in i} (\sum\_{j \in k} e\_{ij})^2} / n_i\\ over the
provider's patients \\k\\, with \\e\_{ij}\\ the residuals: the sandwich
variance of \\\bar Y_i - \bar Z_i^\top\hat\beta\\.

## Examples

``` r
d <- sim_kecc(n_providers = 20, type = "continuous", seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_linear_fe(pd)
test(fit, variance = "robust")
#> Provider tests: wald statistic, robust variance; null median = -0.06783; two.sided, level 0.95
#>   20 providers: 7 above, 8 below, 5 not different, 0 not tested
#>  provider estimate      se null_value statistic         p flag   lower   upper
#>      P001  -0.5629 0.16557   -0.06783    -2.990 2.787e-03   -1 -0.8874 -0.2384
#>      P002  -0.3646 0.08047   -0.06783    -3.688 2.262e-04   -1 -0.5223 -0.2069
#>      P003   0.3252 0.24915   -0.06783     1.578 1.146e-01    0 -0.1631  0.8136
#>      P004   0.4046 0.11173   -0.06783     4.229 2.352e-05    1  0.1856  0.6236
#>      P005   0.4156 0.12481   -0.06783     3.873 1.074e-04    1  0.1710  0.6602
#>      P006  -0.3679 0.08349   -0.06783    -3.594 3.252e-04   -1 -0.5316 -0.2043
#>   ... 14 more rows
```
