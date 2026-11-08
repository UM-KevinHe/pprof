# pprof

<!-- badges: start -->
[![R-CMD-check](https://github.com/UM-KevinHe/pprof/actions/workflows/R-CMD-check.yaml/badge.svg?branch=Yubo)](https://github.com/UM-KevinHe/pprof/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**pprof** profiles health care providers: it fits risk-adjustment models with
one effect per provider, turns the fit into standardized measures, tests every
provider against a reference, calibrates the tests with an empirical null,
measures reliability, selects risk factors, and draws the displays a measure
report needs.

The package is organised by **analysis step**, not by quality measure. A
measure, whatever the outcome it is built on, is assembled from the same
components:

```
pp_data()  →  fit_*()  →  standardize()  →  test()  →  null_model() / calibrate()  →  plot()
  data         model       measures          flags       reference distribution       displays
```

Every step is one function with S3 methods for the model classes, so the same
`test()` call works for a logistic, a linear or a Cox fit, and the documentation
of a step applies to every outcome type.

## Installation

```r
# install.packages("remotes")
remotes::install_github("UM-KevinHe/pprof", ref = "Yubo")
```

A C++17 compiler is required (Rtools on Windows, Xcode command-line tools on
macOS). OpenMP is used when available.

## A complete analysis in ten lines

Patients nested in providers, a binary outcome per patient-month, a few risk
factors:

```r
library(pprof)

d  <- sim_kecc(n_providers = 200, seed = 1)           # or your own data frame
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"),     # declare the columns
              provider = "provider", patient = "patient", cutoff = 10)

fit <- fit_logistic_fe(pd)                            # logistic model, one effect per provider
summary(fit)                                          # the risk-factor coefficients

sm  <- standardize(fit, type = "indirect")            # observed / expected per provider
tst <- test(fit, statistic = "score", variance = "robust",   # one test per provider,
            null_model = "empirical")                 # calibrated by an empirical null
tst                                                   # flags: 1 above, -1 below, 0 not different

plot(tst)                                             # caterpillar plot
plot(tst, type = "funnel")                            # funnel plot with the test's limits
plot(standardize(fit), test = tst)                    # observed against expected
```

## Which functions for which outcome

| Outcome | Example measures | Model | Standardized measure | Tests of the provider effects |
|---|---|---|---|---|
| binary (0/1 per row) | readmission within 30 days of a discharge; a fistula in use in a patient-month | `fit_logistic_fe()` (SerBIN, Firth), `fit_logistic_re()`, `fit_logistic_cre()` | indirect ratio or rate; direct rate | exact, score, Wald × naive, robust, stabilized; likelihood ratio |
| binomial (events out of *n* trials per row) | hospital days out of days at risk in a patient-month | the logistic fitters with `n` | as binary | as binary |
| continuous | a laboratory value, a length of stay | `fit_linear_fe()`, `fit_linear_re()`, `fit_linear_cre()` | indirect and direct differences | score, Wald × naive, robust |
| time to event (counting-process rows) | death; hospital admissions as recurrent events | `fit_cox()` (stratified by provider, two-stage) | observed / expected events, indirect and direct | mid-p Poisson, exact Poisson, log-ratio Wald |
| any of the above, choosing risk factors | measure development | `var_select()` (lasso, group lasso), `cv()` | | |

Every test result can be calibrated with an empirical null (`null_model()`,
`calibrate()`), and every model's `standardize()` result can be tested at the
measure level with `test_measure()`. Reliability is `iur()`.

## Documentation

The site <https://um-kevinhe.github.io/pprof/> has a *Get started* page with
the core concepts, one article per analysis step (data, risk-adjustment models,
variable selection, standardized measures, tests, reference distributions,
reliability, plots), and the function reference grouped the same way.

## Background

The methods are those of the KECC dialysis quality measures and the papers
behind them: block-wise Newton estimation of thousands of provider effects
(He et al. 2013; Wu et al. 2022), Firth-corrected estimation for sparse
providers, patient-clustered and stabilized variances (Pan 2001), exact
Poisson-binomial tests, the empirical null (Efron 2004; Kalbfleisch and Wolfe
2013), direct and indirect standardization (He 2019), inter-unit reliability
(He et al. 2019), penalized selection with provider effects kept, and the
two-stage stratified Cox model of the survival measures (He and Schaubel 2014).
See `citation("pprof")` and the references in each help page.

This is version 2.0.0 (development), which supersedes `pprof` 1.0.3, `pprofR`
and `grplasso`; `NEWS.md` records what changed and why.
