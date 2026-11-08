# Cross-validation of a selection path

Chooses the penalty of a
[`var_select()`](https://um-kevinhe.github.io/pprof/reference/var_select.md)
path by K-fold cross-validation of the binomial deviance (for a survival
outcome, the partial-likelihood deviance per event of the held-out rows
within their own strata and risk sets for the Cox model, the binomial
deviance per person-period for the discrete-time model). Folds are
formed within every provider over its patients (rows when there is no
patient column), so that each provider is present in every training set
whenever it has at least two units; a provider held out entirely is
predicted with the median provider effect of the training fit.

## Usage

``` r
cv(object, ...)

# S3 method for class 'pp_lasso'
cv(object, n_folds = 10, folds = NULL, seed = NULL, ...)
```

## Arguments

- object:

  a path from
  [`var_select()`](https://um-kevinhe.github.io/pprof/reference/var_select.md);
  every fold refits it with the same settings and penalties.

- ...:

  unused.

- n_folds:

  number of folds. Default 10.

- folds:

  optional fold membership (integers from 1) of the rows of `data$data`;
  rows of excluded providers are ignored.

- seed:

  optional seed for the fold assignment; the random-number state of the
  session is restored afterwards.

## Value

An object of class `pp_cv`:

- lambda:

  the penalty sequence of the full-data path.

- cve, cvse:

  the cross-validated deviance per row, event or person-period (mean
  over folds of the fold means) and its standard error across folds.

- lambda_min, lambda_1se:

  the penalty with the smallest `cve` and the largest penalty whose
  `cve` is within one standard error of it.

- selected_min, selected_1se:

  the covariates selected at those penalties.

- fit:

  the path on the full data.

- folds:

  the fold of each row of the fitted data.

## Examples

``` r
# \donttest{
d <- sim_kecc(n_providers = 20, n_covariates = 6, seed = 1)
pd <- pp_data(d, y = "y", x = paste0("x", 1:6), provider = "provider",
              patient = "patient", message = FALSE)
cvp <- cv(var_select(pd, n_lambda = 20), n_folds = 5, seed = 1)
cvp
#> Cross-validated lasso path (5 folds)
#>   lambda_min 0.001232: 6 covariates selected (x1, x2, x3, x4, x5, x6)
#>   lambda_1se 0.02258: 6 covariates selected
cvp$selected_min
#> [1] "x1" "x2" "x3" "x4" "x5" "x6"
# }
```
