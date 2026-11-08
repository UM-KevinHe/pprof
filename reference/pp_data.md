# Prepare data for a provider model

Checks and arranges the data every model and test in the package reads:
the outcome, the covariates, the provider identifier and, when the data
have repeated observations, the patient identifier and the number of
trials; for a survival outcome the time, event and entry columns. Rows
are sorted by provider and patient, provider sizes are counted, and
providers below the size `cutoff` or with no or all events are flagged.

## Usage

``` r
pp_data(
  data,
  y = NULL,
  x,
  provider,
  patient = NULL,
  n = NULL,
  time = NULL,
  event = NULL,
  start = NULL,
  weights = NULL,
  offset = NULL,
  outcome = c("auto", "binary", "binomial", "continuous", "survival"),
  cutoff = 10,
  check = TRUE,
  message = TRUE
)
```

## Arguments

- data:

  a data frame.

- y:

  name of the outcome column: 0/1 for a binary outcome, the number of
  events out of `n` trials for a binomial outcome, or any numeric value
  for a continuous outcome. `NULL` for a survival outcome, which is
  given by `time` and `event`.

- x:

  names of the covariate columns. Covariates must be numeric; code
  factors as 0/1 columns first.

- provider:

  name of the provider identifier column.

- patient:

  name of the patient identifier column, or `NULL` when every row is an
  independent unit. Required by the patient-clustered tests.

- n:

  name of the column with the number of trials per row (binomial
  outcome), or `NULL`.

- time, event:

  for a survival outcome, the names of the columns with the end of each
  row's at-risk interval and its 0/1 event indicator.

- start:

  optional name of the column with the start of the at-risk interval
  (left truncation or counting-process rows); 0 when `NULL`.

- weights:

  optional name of a column of positive case weights, used by
  [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md).

- offset:

  optional name of a column with a known term of the linear predictor,
  used by
  [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md).

- outcome:

  `"auto"` (default), `"binary"`, `"binomial"`, `"continuous"` or
  `"survival"`. With `"auto"` the outcome is survival when `time` or
  `event` is given, binomial when `n` is given, binary when every value
  of `y` is 0 or 1, and continuous otherwise.

- cutoff:

  providers with fewer than `cutoff` patients (rows when `patient` is
  `NULL`) are kept in the data but marked `included = FALSE` and left
  out of model fitting. Default 10. Use 0 to fit every provider, as the
  two-stage survival models of the KECC measures do, and select the
  providers to report at the test stage.

- check:

  run
  [`pp_check()`](https://um-kevinhe.github.io/pprof/reference/pp_check.md)
  first (missing values and zero-variance covariates are errors;
  near-zero variance, high correlation and large VIF are warnings).
  Default `TRUE`.

- message:

  print a short summary. Default `TRUE`.

## Value

An object of class `pp_data`: a list with

- data:

  the sorted data frame with the named columns plus the logical columns
  `included`, `no_events` and `all_events` (the last two are `FALSE` for
  a continuous outcome, and `all_events` for a survival outcome).

- vars:

  the column names `y`, `x`, `provider`, `patient`, `n`, `time`,
  `start`, `event`, `weights`, `offset` and the `outcome` type.

- sizes:

  one row per provider in the sorted order: `provider`, `rows`,
  `patients` (`NA` without a patient column), `included`, `no_events`,
  `all_events`.

- cutoff:

  the cutoff used.

## Details

Rows are ordered by `provider`, then `patient`, with
[`order()`](https://rdrr.io/r/base/order.html)'s radix method, which
sorts character identifiers in the C locale and is stable: a patient's
rows keep their input order, so sort the data by time within patient
first when the stabilized tests will be used. Provider identifiers are
converted to character.

Only the named columns are kept.

## See also

[`pp_check()`](https://um-kevinhe.github.io/pprof/reference/pp_check.md),
[`sim_kecc()`](https://um-kevinhe.github.io/pprof/reference/sim_kecc.md)

## Examples

``` r
d <- sim_kecc(n_providers = 20, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", cutoff = 10, message = FALSE)
pd
#> Prepared data: 2679 rows, 20 providers, 759 patients
#>   outcome y (binary); 3 covariates
#>   20 providers with at least 10 patients included; 0 excluded
#>   included providers with no events: 0; with all events: 0
#>   event rate among included rows: 25.2%
head(pd$sizes)
#>   provider rows patients included no_events all_events
#> 1     P001   64       18     TRUE     FALSE      FALSE
#> 2     P002  197       53     TRUE     FALSE      FALSE
#> 3     P003   45       15     TRUE     FALSE      FALSE
#> 4     P004  164       48     TRUE     FALSE      FALSE
#> 5     P005  122       37     TRUE     FALSE      FALSE
#> 6     P006  208       57     TRUE     FALSE      FALSE

s <- sim_kecc(n_providers = 20, type = "survival", seed = 1)
ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", time = "stop", start = "start",
              event = "event", cutoff = 0, message = FALSE)
ps
#> Prepared data: 2679 rows, 20 providers, 759 patients
#>   survival outcome: time stop from start, event event; 3 covariates
#>   20 providers with at least 0 patients included; 0 excluded
#>   included providers with no events: 0
#>   events among included rows: 336 (12.5% of rows)
```
