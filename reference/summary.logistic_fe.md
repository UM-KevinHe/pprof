# Tests of the covariate coefficients of a logistic fixed-effect fit

Tests of the covariate coefficients of a logistic fixed-effect fit

## Usage

``` r
# S3 method for class 'logistic_fe'
summary(
  object,
  test = c("wald", "robust_wald", "lr", "score"),
  covariates = NULL,
  level = 0.95,
  null = 0,
  ...
)
```

## Arguments

- object:

  a model fitted by
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md).

- test:

  `"wald"` (default; model-based), `"robust_wald"` (patient-clustered
  sandwich; needs a patient column), `"lr"` (likelihood ratio; refits
  the model without the covariate) or `"score"` (refits the model
  without the covariate).

- covariates:

  optional subset of covariate names.

- level:

  confidence level for the Wald intervals.

- null:

  the null value of each coefficient; `"lr"` and `"score"` only allow 0.

- ...:

  unused.

## Value

A data frame with one row per covariate: `covariate`, `estimate`, `se`,
`statistic` (a z statistic for the Wald tests, a chi-square statistic
with one degree of freedom for `"lr"` and `"score"`), `p`, `lower` and
`upper` (Wald tests only).

## Details

The refits use the settings of `object`, including the Firth penalty;
with the penalty the likelihood ratio compares the penalized
log-likelihoods.

## Examples

``` r
d <- sim_kecc(n_providers = 20, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
object <- fit_logistic_fe(pd)
summary(object, test = "robust_wald")
#>   covariate   estimate         se statistic            p      lower      upper
#> 1        x1  0.5942121 0.05031266 11.810389 0.000000e+00  0.4956011  0.6928231
#> 2        x2 -0.4515167 0.10631024 -4.247160 2.164969e-05 -0.6598809 -0.2431524
#> 3        x3  0.3261670 0.05133772  6.353361 2.106610e-10  0.2255469  0.4267871
```
