# Simulate data shaped like a KECC quality measure

Patients nested in providers with repeated observations per patient, a
provider effect, a patient random effect and a few covariates. The
outcome is binary per row (a patient-month, say), binomial out of a
number of trials per row (days at risk, say), continuous, or a survival
time with censoring, split into counting-process rows.

## Usage

``` r
sim_kecc(
  n_providers = 50,
  provider_size = c(15, 60),
  rows_per_patient = c(1, 6),
  type = c("binary", "binomial", "continuous", "survival"),
  trials = c(20, 31),
  n_covariates = 3,
  beta = NULL,
  gamma_sd = 0.5,
  patient_sd = 0.5,
  base_rate = 0.25,
  sigma = 1,
  base_hazard = 0.5,
  follow_up = c(0.5, 2),
  seed = NULL
)
```

## Arguments

- n_providers:

  number of providers.

- provider_size:

  range of the number of patients per provider; sizes are drawn
  uniformly between the two values.

- rows_per_patient:

  range of the number of rows per patient.

- type:

  `"binary"` (one trial per row), `"binomial"` (`n` trials per row),
  `"continuous"` or `"survival"`.

- trials:

  range of the number of trials per row for `type = "binomial"`.

- n_covariates:

  number of covariates (at least 1): `x1` is a patient-level standard
  normal, `x2` a patient-level 0/1 with probability 0.4, and any further
  covariates are row-level standard normals (patient-level for a
  survival outcome).

- beta:

  covariate coefficients; by default `0.5, -0.4, 0.3` recycled.

- gamma_sd:

  standard deviation of the provider effects around the logit of
  `base_rate` (binary, binomial) or around 0 (continuous).

- patient_sd:

  standard deviation of the patient random effect.

- base_rate:

  event probability of an average patient at the median provider (binary
  and binomial outcomes).

- sigma:

  residual standard deviation of a continuous outcome.

- base_hazard, follow_up:

  for `type = "survival"`: the hazard of an average patient at the
  median provider (exponential event times), and the range of the
  censoring time, drawn uniformly. Each patient's follow-up is cut into
  `rows_per_patient` consecutive rows `(start, stop]`, the last carrying
  the event if it came before censoring.

- seed:

  optional random seed.

## Value

A data frame with columns `provider`, `patient`, `n` (trials; 1 for
binary and continuous outcomes), `y` and the covariates, one row per
observation, in random order; for a survival outcome the columns
`start`, `stop` and `event` replace `n` and `y`. The true parameters are
kept in `attr(, "truth")`: `gamma` (named by provider), `beta`,
`patient_effect` (named by patient).

## Examples

``` r
d <- sim_kecc(n_providers = 5, seed = 1)
str(attr(d, "truth"))
#> List of 3
#>  $ gamma         : Named num [1:5] -0.572 -1.658 -0.931 -0.851 -1.03
#>   ..- attr(*, "names")= chr [1:5] "P001" "P002" "P003" "P004" ...
#>  $ beta          : Named num [1:3] 0.5 -0.4 0.3
#>   ..- attr(*, "names")= chr [1:3] "x1" "x2" "x3"
#>  $ patient_effect: Named num [1:171] -0.0594 0.0988 -0.5343 -0.4016 -0.5569 ...
#>   ..- attr(*, "names")= chr [1:171] "P001-001" "P001-002" "P001-003" "P001-004" ...
s <- sim_kecc(n_providers = 5, type = "survival", seed = 1)
head(s)
#>   provider  patient     start       stop event         x1 x2         x3
#> 1     P003 P003-010 0.0000000 0.13188650     1  0.1063206  1 -0.2464703
#> 2     P004 P004-027 0.4887299 0.57594207     1  0.4596217  0  0.7048330
#> 3     P002 P002-025 0.8418941 1.06252226     0  0.2222181  0 -0.7610403
#> 4     P002 P002-012 0.5227791 1.03538921     0  0.4345367  1  0.5756203
#> 5     P004 P004-004 0.0000000 0.04786032     0 -0.1296354  0  0.7159321
#> 6     P005 P005-028 0.6503206 0.79878563     0 -0.8846849  1  1.3609245
```
