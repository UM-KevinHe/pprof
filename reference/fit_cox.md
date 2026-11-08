# Cox model with provider strata or a pooled baseline

Fits the Cox proportional hazards model of the survival measures (SMR,
SHR, STrR) through
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html). By
default the provider is a stratification variable, so that each provider
has its own baseline hazard and the covariate effects are estimated
within providers (He and Schaubel 2014); the national baseline is then
the Breslow estimator over every row with the fitted linear predictor as
an offset. This is the two-stage model of the KECC measures: stage one
estimates the covariate effects in the stratified model, stage two
estimates the optional covariates `x2` with the stage-one linear
predictor as offset and gives the baseline. With `strata = "none"` one
pooled model is fitted, without provider effects.

## Usage

``` r
fit_cox(
  data,
  strata = c("provider", "none"),
  x2 = NULL,
  ties = c("breslow", "efron"),
  message = FALSE
)
```

## Arguments

- data:

  an object from
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  with a survival outcome (`time`, `event`, optional `start`, `weights`
  and `offset`).

- strata:

  `"provider"` (default) or `"none"`.

- x2:

  optional names of covariate columns fitted in the second stage with
  the first-stage linear predictor as offset (calendar-year indicators
  in the KECC measures, say); they must be in `data`'s covariates `x` or
  be further numeric columns of the original data frame named in `x`.
  `NULL` (default) for no second-stage covariates.

- ties:

  `"breslow"` (default, the KECC convention) or `"efron"`.

- message:

  print the fit. Default `FALSE`.

## Value

An object of class `cox`: a list with

- coefficients:

  `beta` (first stage) and `beta2` (second stage, or `NULL`).

- stage1, stage2:

  the [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html)
  fits (`stage2` is `NULL` without `x2`).

- linear_predictor:

  the full linear predictor of every included row, including the offset.

- expected, residuals:

  expected events and martingale residuals per row.

- baseline:

  a data frame with the distinct event times and the national cumulative
  baseline hazard.

- data, vars, sizes, rows, control, converged:

  as in
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md).

## Details

Rows with `included = FALSE` are left out; set `cutoff = 0` in
[`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
to fit every provider, as the measures do, and choose the providers to
report in
[`test()`](https://um-kevinhe.github.io/pprof/reference/test.md).

When `data` has a patient column it is passed as `cluster` to
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html), so
that the robust (patient-clustered) covariance is available from
[`vcov.cox()`](https://um-kevinhe.github.io/pprof/reference/vcov.cox.md)
and [`summary()`](https://rdrr.io/r/base/summary.html).

The expected number of events of every row is \\E_i =
\exp(\eta_i)\\\[\Lambda_0(t_i) - \Lambda_0(s_i)\]\\ with \\\eta_i\\ the
full linear predictor (covariates of both stages and the offset) and
\\\Lambda_0\\ the Breslow estimator of the national cumulative baseline
hazard at \\\eta = 0\\, weighted by the case weights; the row's
martingale residual is \\\delta_i - E_i\\. For `ties = "breslow"` these
are the residuals
[`survival::residuals.coxph()`](https://rdrr.io/pkg/survival/man/residuals.coxph.html)
returns for the second-stage model, and `expected` is what the KECC
measures compute as `event - residual`.

## References

He K, Schaubel DE (2014). Standardized mortality ratio for evaluating
center-specific mortality: assessment and alternative. *Statistics in
Biosciences*, 6(2), 1-19.

Liu D, Schaubel DE, Kalbfleisch JD (2012). Computationally efficient
marginal models for clustered recurrent event data. *Biometrics*, 68(2),
637-647.

## See also

[`standardize.cox()`](https://um-kevinhe.github.io/pprof/reference/standardize.cox.md),
[`test.cox()`](https://um-kevinhe.github.io/pprof/reference/test.cox.md),
[`summary.cox()`](https://um-kevinhe.github.io/pprof/reference/summary.cox.md)

## Examples

``` r
s <- sim_kecc(n_providers = 20, type = "survival", seed = 1)
ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", time = "stop", start = "start",
              event = "event", cutoff = 0, message = FALSE)
fit <- fit_cox(ps)
fit
#> Cox model stratified by provider (two-stage)
#>   2679 rows, 20 providers, 3 covariates; 336 events; ties breslow
#>   expected events sum to 336
coef(fit)$beta
#>         x1         x2         x3 
#>  0.5168166 -0.1163487  0.2665547 
```
