# pprof: provider profiling by components

Components for profiling health care providers, organised by analysis
step: data preparation
([`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md),
[`pp_check()`](https://um-kevinhe.github.io/pprof/reference/pp_check.md)),
risk-adjustment models, standardized measures
([`standardize()`](https://um-kevinhe.github.io/pprof/reference/standardize.md)),
tests of provider effects
([`test()`](https://um-kevinhe.github.io/pprof/reference/test.md)) and
covariates ([`summary()`](https://rdrr.io/r/base/summary.html)),
reference distributions, reliability, variable selection and plots.
Every quality measure is assembled from these components; nothing in the
package is specific to one measure.

## Conventions

- Every function takes a data frame and column names (`y`, `x`,
  `provider`, `patient`, `n`; `time`, `start`, `event` for survival
  outcomes).

- Provider identifiers are handled as character strings and rows are
  sorted by provider, then patient, in the C locale; a patient's rows
  keep their input order, which the stabilized tests rely on.

- Rates are proportions on the 0-1 scale.

- Flags are `1` (above the reference), `-1` (below), `0` (not different)
  or `NA` (not tested); which tail is "worse" is a property of the
  measure and is left to the caller.

## See also

Useful links:

- <https://um-kevinhe.github.io/pprof>

- <https://github.com/UM-KevinHe/pprof>

- Report bugs at <https://github.com/UM-KevinHe/pprof/issues>

## Author

**Maintainer**: Yubo Shao <ybshao@umich.edu>

Authors:

- Kevin He <kevinhe@umich.edu>
