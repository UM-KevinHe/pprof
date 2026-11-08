# Fit the linear model with fixed provider effects

Fits \\Y\_{ij} = \gamma_i + Z\_{ij}^\top\beta + \epsilon\_{ij}\\ by
least squares on the providers that
[`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
marked as included. The covariate coefficients come from the
within-provider regression (every variable centred at its provider
mean), so no provider dummy variables are formed; the provider effects
are \\\hat\gamma_i = \bar Y_i - \bar Z_i^\top\hat\beta\\.

## Usage

``` r
fit_linear_fe(data)
```

## Arguments

- data:

  an object from
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  with a continuous (or any numeric) outcome.

## Value

An object of class `linear_fe`:

- coefficients:

  a list with `gamma` (provider effects, named by provider) and `beta`
  (covariate coefficients, named by covariate).

- sigma:

  the residual standard deviation, \\\hat\sigma^2 = RSS / (n - m - p)\\.

- fitted, linear_predictor, residuals:

  per included row; `linear_predictor` is \\Z\hat\beta\\.

- loglik, aic, bic:

  Gaussian log-likelihood at the maximum-likelihood variance and the
  information criteria with \\m + p + 1\\ parameters.

- df_residual:

  \\n - m - p\\.

- vcov_beta:

  \\\hat\sigma^2 (\tilde Z^\top \tilde Z)^{-1}\\ with \\\tilde Z\\ the
  within-centred covariates.

- data, vars, sizes, cutoff, rows:

  as in
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md).

## See also

[`test.linear_fe()`](https://um-kevinhe.github.io/pprof/reference/test.linear_fe.md),
[`summary.linear_fe()`](https://um-kevinhe.github.io/pprof/reference/summary.linear_fe.md),
[`standardize.linear_fe()`](https://um-kevinhe.github.io/pprof/reference/standardize.linear_fe.md),
[`confint.linear_fe()`](https://um-kevinhe.github.io/pprof/reference/confint.linear_fe.md),
[`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md)

## Examples

``` r
d <- sim_kecc(n_providers = 20, type = "continuous", seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_linear_fe(pd)
fit
#> Linear model with fixed provider effects
#>   2679 rows, 20 providers, 3 covariates
#>   residual standard deviation 1.119 on 2656 degrees of freedom
#>   provider effects: median -0.06783, range -0.6820 to  0.8072
```
