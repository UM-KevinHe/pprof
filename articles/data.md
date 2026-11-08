# Data preparation

Every component of the package reads one prepared object, made once by
[`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
from a data frame and the names of its columns. This page explains what
the object records and the decisions it embodies.

## The unit of observation

A row of the data is one **observation**: whatever the measure counts
outcomes on.

| Unit | Example | Outcome per row |
|----|----|----|
| patient-month | one row per patient per calendar month | 0/1 (an event in the month), or a count of days out of the days at risk |
| index discharge | one row per hospital discharge | 0/1 (readmitted within 30 days) |
| patient-record | one row per patient per provider per year | a time at risk and an event indicator |
| patient | one row per patient | a continuous value |

A patient usually contributes several rows, and those rows are **not
independent**: a patient sicker than the risk factors say is sicker in
every month.
[`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
therefore takes the patient identifier as well as the provider
identifier. The patient column is what the patient-clustered (robust and
stabilized) tests need; without it only the naive tests are available.

## Declaring the columns

``` r

library(pprof)
d <- sim_kecc(n_providers = 30, seed = 1)
str(d)
#> 'data.frame':    4343 obs. of  7 variables:
#>  $ provider: chr  "P026" "P027" "P020" "P013" ...
#>  $ patient : chr  "P026-013" "P027-027" "P020-029" "P013-017" ...
#>  $ n       : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ y       : int  0 0 1 0 0 1 1 0 1 0 ...
#>  $ x1      : num  0.2617 -0.0742 -1.5689 1.1322 -1.3084 ...
#>  $ x2      : num  1 0 0 1 1 0 0 0 0 0 ...
#>  $ x3      : num  2.2451 -1.5913 0.0299 -1.1788 -0.8343 ...
#>  - attr(*, "truth")=List of 3
#>   ..$ gamma         : Named num [1:30] 0.429 -1.229 -1.326 -1.02 -0.632 ...
#>   .. ..- attr(*, "names")= chr [1:30] "P001" "P002" "P003" "P004" ...
#>   ..$ beta          : Named num [1:3] 0.5 -0.4 0.3
#>   .. ..- attr(*, "names")= chr [1:3] "x1" "x2" "x3"
#>   ..$ patient_effect: Named num [1:1237] -0.3115 -0.2198 -0.2526 0.093 0.0882 ...
#>   .. ..- attr(*, "names")= chr [1:1237] "P001-001" "P001-002" "P001-003" "P001-004" ...
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", cutoff = 10)
#> Checking missing values ...
#> Checking variation in covariates ...
#> Checking pairwise correlation among covariates ...
#> Checking variance inflation factors ...
#> Checks completed.
#> Prepared data: 4343 rows, 30 providers, 1237 patients
#>   outcome y (binary); 3 covariates
#>   30 providers with at least 10 patients included; 0 excluded
#>   included providers with no events: 0; with all events: 0
#>   event rate among included rows: 25.1%
```

- `y`: the outcome. Binary (0/1), binomial (the number of events out of
  the `n` trials of the row, for instance hospital days out of days at
  risk) or continuous; `outcome = "auto"` tells them apart, and
  `outcome` can be given explicitly.
- `x`: the risk factors, numeric. Code a factor as 0/1 columns first
  (the package does not expand factors, so that the design matrix is
  exactly what you declared).
- `provider`, `patient`: identifiers, converted to character.
- `n`: the trials column of a binomial outcome.

### Survival outcomes

A time-to-event outcome is given by its at-risk interval and event
indicator instead of `y`:

``` r

s <- sim_kecc(n_providers = 30, type = "survival", seed = 1)
head(s)
#>   provider  patient     start      stop event            x1 x2           x3
#> 1     P024 P024-040 0.0000000 0.5346169     0 -1.3394024080  0 -2.758961290
#> 2     P014 P014-023 0.2534563 0.5644242     0  0.0711682667  0 -0.550372599
#> 3     P004 P004-025 0.5315569 1.5701286     0  2.1736849196  0  0.002432616
#> 4     P002 P002-034 1.7584042 1.7838728     0 -0.4243772091  0 -1.640790694
#> 5     P020 P020-006 0.0000000 0.2961319     1  0.0005632156  1  0.364351557
#> 6     P009 P009-032 0.0000000 0.1432058     0 -0.4120366654  1 -0.939521327
ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider", patient = "patient",
              time = "stop", start = "start", event = "event", cutoff = 0)
#> Checking missing values ...
#> Checking variation in covariates ...
#> Checking pairwise correlation among covariates ...
#> Checking variance inflation factors ...
#> Checks completed.
#> Prepared data: 4343 rows, 30 providers, 1237 patients
#>   survival outcome: time stop from start, event event; 3 covariates
#>   30 providers with at least 0 patients included; 0 excluded
#>   included providers with no events: 0
#>   events among included rows: 549 (12.6% of rows)
```

`time` is the end of the interval, `start` its beginning (0 when absent:
no left truncation), `event` 0/1. A patient’s follow-up may be cut into
several rows (counting-process form), for instance one per calendar year
or one per provider the patient was attributed to. Two more columns are
read by the Cox model only: `weights` (case weights) and `offset` (a
known term of the linear predictor).

## What `pp_data()` records

``` r

pd
#> Prepared data: 4343 rows, 30 providers, 1237 patients
#>   outcome y (binary); 3 covariates
#>   30 providers with at least 10 patients included; 0 excluded
#>   included providers with no events: 0; with all events: 0
#>   event rate among included rows: 25.1%
head(pd$sizes)
#>   provider rows patients included no_events all_events
#> 1     P001   62       18     TRUE     FALSE      FALSE
#> 2     P002  193       53     TRUE     FALSE      FALSE
#> 3     P003   51       15     TRUE     FALSE      FALSE
#> 4     P004  158       48     TRUE     FALSE      FALSE
#> 5     P005  129       37     TRUE     FALSE      FALSE
#> 6     P006  214       57     TRUE     FALSE      FALSE
pd$vars
#> $y
#> [1] "y"
#> 
#> $x
#> [1] "x1" "x2" "x3"
#> 
#> $provider
#> [1] "provider"
#> 
#> $patient
#> [1] "patient"
#> 
#> $n
#> NULL
#> 
#> $time
#> NULL
#> 
#> $start
#> NULL
#> 
#> $event
#> NULL
#> 
#> $weights
#> NULL
#> 
#> $offset
#> NULL
#> 
#> $outcome
#> [1] "binary"
```

- `data`: the named columns, sorted by provider and then patient. The
  sort is stable and in the C locale, so a patient’s rows keep their
  input order; sort by time within patient beforehand when the
  stabilized variance will be used.
- `sizes`: one row per provider with its number of rows and patients,
  and three flags: `included`, `no_events`, `all_events`.
- `vars`: the column names and the outcome type, which every later
  function reads, so that nothing is declared twice.

## The cutoff

Providers with fewer than `cutoff` patients (rows, when there is no
patient column) are kept in the data with `included = FALSE` and left
out of model fitting: too little information to estimate an effect, and
in the measures such providers are reported as “too few patients”. The
comparison is `>=` on the number of patients.

For the two-stage survival models the convention is different: every
provider enters the stratified fit (`cutoff = 0`) and the providers to
report are chosen at the test stage (`test(fit, min_expected = 3)` and
the `providers` argument).

## Providers whose outcomes are all the same

A provider with no events, or with an event in every row, has no finite
maximum-likelihood effect: `no_events` and `all_events` mark them. Such
providers are fitted at the bound of the parameter space (`bound` in
[`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md),
10 on the logit scale around the median), they can be tested with the
score or exact statistic but not with the Wald statistic, and the Firth
penalty (`fit_logistic_fe(penalty = "firth")`) gives them a finite
estimate.

## Checks

[`pp_check()`](https://um-kevinhe.github.io/pprof/reference/pp_check.md)
runs first (`check = TRUE`): missing values in the named columns and
constant risk factors are errors; near-constant risk factors, pairs of
risk factors with correlation above 0.9 and variance inflation factors
above 10 are warnings. It can be called on its own with other
thresholds:

``` r

pp_check(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider", cor_threshold = 0.8)
#> Checking missing values ...
#> Checking variation in covariates ...
#> Checking pairwise correlation among covariates ...
#> Checking variance inflation factors ...
#> Checks completed.
```

## Simulated data

[`sim_kecc()`](https://um-kevinhe.github.io/pprof/reference/sim_kecc.md)
simulates data of every outcome type with the structure of the measures:
patients nested in providers, several rows per patient, a provider
effect, a patient random effect that makes a patient’s rows correlated,
and a few risk factors. The true parameters are kept in
`attr(d, "truth")`. The articles use it throughout; replace it with your
data frame.
