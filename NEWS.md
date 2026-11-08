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

## Still to come (in this order)

* The logistic fixed-effect engine of `pprofR` under the unified interface.
* Random, correlated-random and linear models of 1.0.3 with their defects
  fixed.
* Variable selection of `grplasso` with its defects fixed.
* Empirical null, inter-unit reliability, Cox O/E components.
* Plots.
