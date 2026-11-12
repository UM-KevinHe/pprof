# Fit the discrete-time hazard model with fixed provider effects

Fits \\\mathrm{logit}\\ h\_{ijk} = \alpha_k + \gamma_i +
Z\_{ij}^\top\beta\\, the hazard of the event in period \\k\\ for subject
\\j\\ of provider \\i\\, by maximum likelihood on one row per subject.
The person-periods are looped over inside the computation and never
stored.

## Usage

``` r
fit_discrete_fe(
  data,
  max_iter = 10000,
  tol = 1e-05,
  bound = 10,
  backtrack = TRUE,
  message = FALSE
)
```

## Arguments

- data:

  an object from
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  with a survival outcome: `time` is the subject's last period (the
  event period, or the last period observed without the event) and
  `event` the event indicator. One row per subject, without `start`,
  `weights` or `offset`.

- max_iter:

  maximum number of Newton iterations. Default 10,000.

- tol:

  convergence tolerance on the period effects and covariate
  coefficients. Default 1e-5.

- bound:

  provider effects are kept within `median(gamma) +/- bound`. Default
  10.

- backtrack:

  use a backtracking line search. Default `TRUE`.

- message:

  print the iteration history.

## Value

An object of class `discrete_fe`:

- coefficients:

  a list with `gamma` (provider effects), `alpha` (period effects, named
  by period, the first fixed at zero) and `beta` (covariate
  coefficients).

- linear_predictor:

  \\Z\hat\beta\\ of the included subjects.

- loglik, aic, bic:

  log-likelihood and information criteria (the BIC counts
  person-periods).

- periods:

  the number of periods each subject is at risk.

- time_grid:

  the event periods.

- person_periods:

  the total number of person-periods at risk.

- data, vars, sizes, cutoff, rows, control, converged:

  as for
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md);
  `rows` counts subjects.

## Details

The periods are the distinct values of `time` among the events: a period
in which no subject has the event has a hazard of zero and adds nothing
to the likelihood. A subject is at risk in every period up to its own
time and has the event, if any, in its last one. The first period's
effect is fixed at zero, so \\\gamma_i\\ is the log-odds of the hazard
in the first period at \\Z = 0\\; tests and standardized measures do not
depend on this choice.

The likelihood is that of the logistic model on the person-period data,
and every result (coefficients, variances, tests, standardized measures)
equals what
[`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md)
gives on the expanded data with the period effects as covariates. The
Newton step solves through the Schur complement of the diagonal provider
block of the information, as SerBIN does.

## See also

[`test.discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/test.discrete_fe.md),
[`standardize.discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/standardize.discrete_fe.md),
[`confint.discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/confint.discrete_fe.md),
[`summary.discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/summary.discrete_fe.md),
and `var_select(model = "discrete")` for selecting the covariates.

## Examples

``` r
s <- sim_kecc(n_providers = 20, provider_size = c(20, 40), rows_per_patient = c(1, 1),
              type = "survival", seed = 1)
s$month <- pmin(ceiling(s$stop * 12), 12)
s$event[s$stop * 12 > 12] <- 0
pd <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider", patient = "patient",
              time = "month", event = "event", message = FALSE)
fit <- fit_discrete_fe(pd)
fit
#> Discrete-time hazard model with fixed provider effects
#>   599 subjects, 5355 person-periods, 20 providers, 12 event periods, 3 covariates
#>   log-likelihood -880.81; iterations 5
#>   provider effects: median -3.013, range -3.907 to -1.969
coef(fit)$alpha
#>            1            2            3            4            5            6 
#>  0.000000000 -0.425938274 -0.001920639  0.003559516  0.213482472 -0.055165573 
#>            7            8            9           10           11           12 
#>  0.065604829 -0.218351045 -0.558239003  0.050250679  0.053689440 -0.619540070 
```
