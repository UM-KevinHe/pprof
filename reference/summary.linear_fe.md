# Tests of the covariate coefficients of a linear fixed-effect fit

Tests of the covariate coefficients of a linear fixed-effect fit

## Usage

``` r
# S3 method for class 'linear_fe'
summary(
  object,
  test = c("wald", "robust_wald"),
  covariates = NULL,
  level = 0.95,
  null = 0,
  ...
)
```

## Arguments

- object:

  a model fitted by
  [`fit_linear_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_linear_fe.md).

- test:

  `"wald"` (default; model-based t test with \\n - m - p\\ degrees of
  freedom) or `"robust_wald"` (patient-clustered sandwich, standard
  normal reference; needs a patient column).

- covariates:

  optional subset of covariate names.

- level:

  confidence level for the intervals.

- null:

  the null value of each coefficient.

- ...:

  unused.

## Value

A data frame with one row per covariate: `covariate`, `estimate`, `se`,
`statistic`, `p`, `lower`, `upper`.

## Examples

``` r
d <- sim_kecc(n_providers = 20, type = "continuous", seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
object <- fit_linear_fe(pd)
summary(object)
#>    covariate   estimate         se statistic             p      lower
#> x1        x1  0.5201115 0.02041176 25.480973 2.563444e-128  0.4800870
#> x2        x2 -0.3800450 0.04434541 -8.570108  0.000000e+00 -0.4670000
#> x3        x3  0.3101395 0.02133061 14.539643  3.838339e-46  0.2683132
#>         upper
#> x1  0.5601361
#> x2 -0.2930899
#> x3  0.3519658
```
