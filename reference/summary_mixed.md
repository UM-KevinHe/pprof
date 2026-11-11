# Tests of the fixed effects of a random-effect fit

Wald tests of the fixed-effect coefficients (including the intercept)
with the covariance matrix of the lme4 object, referred to the standard
normal.

## Usage

``` r
# S3 method for class 'logistic_re'
summary(object, covariates = NULL, level = 0.95, null = 0, ...)

# S3 method for class 'linear_re'
summary(object, covariates = NULL, level = 0.95, null = 0, ...)
```

## Arguments

- object:

  a model from
  [`fit_logistic_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md),
  [`fit_logistic_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md),
  [`fit_linear_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md)
  or
  [`fit_linear_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md).

- covariates:

  optional subset of term names (as in `names(coef(object)$beta)`).

- level:

  confidence level for the intervals.

- null:

  the null value of each coefficient.

- ...:

  unused.

## Value

A data frame with one row per term: `covariate`, `estimate`, `se`,
`statistic`, `p`, `lower`, `upper`.

## Examples

``` r
# \donttest{
d <- sim_kecc(n_providers = 20, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
summary(fit_logistic_re(pd))
#>               covariate   estimate         se statistic            p      lower
#> (Intercept) (Intercept) -1.0042891 0.13013167 -7.717484 1.176836e-14 -1.2593424
#> x1                   x1  0.5839857 0.04814198 12.130486 7.281575e-34  0.4896291
#> x2                   x2 -0.4541446 0.09981072 -4.550059 5.363096e-06 -0.6497700
#> x3                   x3  0.3245055 0.04780162  6.788589 1.132355e-11  0.2308161
#>                  upper
#> (Intercept) -0.7492357
#> x1           0.6783422
#> x2          -0.2585192
#> x3           0.4181950
# }
```
