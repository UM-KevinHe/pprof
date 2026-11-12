# Variance of the provider effects for a Wald test

Two model-based variances of a fixed provider effect \\\hat\gamma_i\\
are available. The `"case_mix"` variance is that of \\\hat\gamma_i +
\bar x_i^\top\hat\beta\\, the provider's outcome (log-odds or mean) at
its own average case mix \\\bar x_i\\; it does not depend on where the
covariates are centred and is the default of
[`test()`](https://um-kevinhe.github.io/pprof/reference/test.md). The
`"full"` variance is the \\(\gamma_i,\gamma_i)\\ element of the inverse
information, \\1/I\_{ii} + \bar z_i^\top \mathrm{Var}(\hat\beta)\bar
z_i\\; it is the variance of the outcome at \\Z = 0\\ and changes with
the origin of the covariates, so it is only appropriate when the null
value is an external constant stated in the same coding.

## Usage

``` r
# S3 method for class 'discrete_fe'
gamma_variance(fit, type = c("case_mix", "full"), ...)

# S3 method for class 'linear_fe'
gamma_variance(fit, type = c("case_mix", "full"), ...)

gamma_variance(fit, type = c("case_mix", "full"), ...)

# S3 method for class 'logistic_fe'
gamma_variance(fit, type = c("case_mix", "full"), ...)
```

## Arguments

- fit:

  a model fitted by
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md)
  or
  [`fit_linear_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_linear_fe.md).

- type:

  `"case_mix"` (default) or `"full"`.

- ...:

  unused.

## Value

A named numeric vector, one variance per provider.

## Details

For the logistic model \\\bar z_i\\ is the \\N p (1 - p)\\-weighted mean
covariate and the case-mix variance is \\1/I\_{ii} + (\bar z_i - \bar
x_i)^\top \mathrm{Var}(\hat\beta) (\bar z_i - \bar x_i)\\; for the
linear model \\\bar z_i = \bar x_i\\ and the case-mix variance is
exactly \\\hat\sigma^2 / n_i\\.
