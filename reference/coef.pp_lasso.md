# Coefficients and predictions along a selection path

Coefficients and predictions along a selection path

## Usage

``` r
# S3 method for class 'pp_lasso'
coef(object, lambda = NULL, which = NULL, drop = TRUE, ...)

# S3 method for class 'pp_lasso'
predict(
  object,
  newdata = NULL,
  lambda = NULL,
  which = NULL,
  type = c("link", "response"),
  unseen = c("missing", "median"),
  ...
)

selected_covariates(object, lambda = NULL, which = NULL)
```

## Arguments

- object:

  a path from
  [`var_select()`](https://um-kevinhe.github.io/pprof/reference/var_select.md)
  or
  [`var_select()`](https://um-kevinhe.github.io/pprof/reference/var_select.md).

- lambda:

  penalty values at which to report; values between two points of the
  path are interpolated linearly.

- which:

  alternatively, indices into `object$lambda`; default all.

- drop:

  return vectors when a single value is requested.

- ...:

  unused.

- newdata:

  a data frame with the covariate and provider columns of the fit;
  default the fitted rows.

- type:

  `"link"` (the linear predictor, including any offset) or `"response"`:
  the probability for the logistic model, the relative hazard \\e^\eta\\
  for the Cox model, and for the discrete-time model the probability of
  the event by the row's time, \\1 - \prod_k (1 - h\_{jk})\\ over the
  fit's periods up to that time (a time between two periods counts up to
  the earlier one, a time before the first period counts the first).

- unseen:

  for providers of `newdata` that are not in the fit: `"missing"`
  (default; `NA`) or `"median"` (the median provider effect). Not used
  by the Cox model, which has no provider parameters.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html): a list with `gamma`
(`NULL` for the Cox model) and `beta`, plus `alpha` for the
discrete-time model (matrices, or vectors with `drop = TRUE` and one
value of `lambda`). [`predict()`](https://rdrr.io/r/stats/predict.html):
a matrix with one column per requested `lambda`.
`selected_covariates()`: the names of the penalized covariates with a
nonzero coefficient, as a list over `lambda`.

## Examples

``` r
d <- sim_kecc(n_providers = 20, n_covariates = 6, seed = 1)
pd <- pp_data(d, y = "y", x = paste0("x", 1:6), provider = "provider",
              patient = "patient", message = FALSE)
path <- var_select(pd, n_lambda = 20)
selected_covariates(path, which = c(5, 15))
#> $`0.0139041`
#> [1] "x1" "x2" "x3" "x4" "x5" "x6"
#> 
#> $`0.000109114`
#> [1] "x1" "x2" "x3" "x4" "x5" "x6"
#> 
head(predict(path, which = 15, type = "response"))
#>      0.000109114
#> [1,]   0.5765182
#> [2,]   0.4244555
#> [3,]   0.5604338
#> [4,]   0.4192023
#> [5,]   0.2504631
#> [6,]   0.2823876
```
