# pprof 2.0.0 (development)

A rewrite of `pprof` 1.0.3 that absorbs `pprofR` 0.2.0 and `grplasso` 1.0.
The package is organised by analysis step, and every function takes a data
frame with column names.

## Data

* `pp_data()` replaces `fe_data_prep()` and the data handling inside the
  1.0.3 fitters: one prepared object (rows sorted by provider and patient in
  the C locale, patient rows kept in input order, provider sizes in patients
  or rows, `cutoff` with `>=`, flags for providers with no or all events).
* `pp_check()` replaces `data_check()`: no `caret` or `olsrr`, thresholds as
  arguments, metrics returned.
* `sim_kecc()` simulates data shaped like the KECC measures (patient-month
  binary or binomial days outcomes, patient random effect).

## Logistic fixed-effect model (the `pprofR` engine under the unified interface)

* `fit_logistic_fe()` replaces `logis_fe()`: binary or binomial outcomes, SerBIN or BAN.
* `test_providers()` replaces `test_fe()`: `statistic` (exact, score, Wald) x `variance` (naive,
  robust, stabilized) with one fixed set of output columns; the naive Wald statistic uses the
  case-mix variance (`gamma_variance()`), which does not depend on the covariate origin; the
  restricted score test of 1.0.3 is back as `score_type = "restricted"`.
* `test_covariates()` replaces `summary_fe_covar()`; the likelihood-ratio refit respects the fit's
  own `bound`.
* `confint()` gives intervals for the provider effects; `standardize()` replaces `SR_output()` and
  maps those limits to the ratios and rates; `direct_expectation()` replaces `direct_exp()`.
* `test_measure()` replaces `test_SM()` and adds the log scale.

## Still to come (in this order)

* Random, correlated-random and linear models of 1.0.3 with their defects
  fixed.
* Variable selection of `grplasso` with its defects fixed.
* Empirical null, inter-unit reliability, Cox O/E components.
* Plots.
