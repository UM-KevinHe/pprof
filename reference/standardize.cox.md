# Standardized ratios of a Cox fit

Indirect or direct standardization of every provider's event count from
a Cox model, in the form of the survival measures (SMR, SHR).

## Usage

``` r
# S3 method for class 'cox'
standardize(fit, type = c("indirect", "direct"), providers = NULL, ...)
```

## Arguments

- fit:

  a model fitted by
  [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md).

- type:

  `"indirect"` (default) or `"direct"`.

- providers:

  optional subset of provider identifiers; the baselines always use
  every fitted row.

- ...:

  unused.

## Value

A data frame with one row per provider. Indirect: `provider`,
`observed`, `expected`, `ratio`, `person_time`. Direct: `provider`,
`observed` (the total \\O\\), `expected` (\\E^{(j)}\\), `ratio`, `n_pop`
(the number of fitted rows). Attribute `type`.

## Details

With \\w_i\\ the case weights, \\\delta_i\\ the event indicator, \\(s_i,
t_i\]\\ the at-risk interval and \\\eta_i\\ the fitted linear predictor
of row \\i\\:

- indirect: \\O_j = \sum\_{i \in j} w_i \delta_i\\ and \\E_j = \sum\_{i
  \in j} w_i\\ e^{\eta_i} \[\Lambda_0(t_i) - \Lambda_0(s_i)\]\\, the
  events expected for provider \\j\\'s rows at the national baseline
  \\\Lambda_0\\ (the Breslow estimator over every row); ratio \\O_j /
  E_j\\; the expected counts add up to the observed ones;

- direct: \\E^{(j)} = \sum_i w_i\\ e^{\eta_i} \[\Lambda\_{0j}(t_i) -
  \Lambda\_{0j}(s_i)\]\\ over every row, the events expected if every
  patient had provider \\j\\'s own baseline \\\Lambda\_{0j}\\ (the
  Breslow estimator over provider \\j\\'s rows at the same
  coefficients); ratio \\E^{(j)} / O\\ with \\O\\ the total observed
  events. Each Breslow increment of provider \\j\\ being one weighted
  event over its risk-set sum \\RS_j(t)\\, \\E^{(j)}\\ equals \\\sum\\
  over provider \\j\\'s events of \\w_i\\ RS(t_i) / RS_j(t_i)\\, the
  ratio of the population's to the provider's risk-set sums of \\w
  e^\eta\\.

Both baselines are Breslow estimators whatever `ties` fitted the
coefficients. The person-time is \\\sum\_{i \in j} (t_i - s_i)\\, the
grouping variable of the survival measures' empirical null.

## References

He K, Schaubel DE (2014). Standardized mortality ratio for evaluating
center-specific mortality: assessment and alternative. *Statistics in
Biosciences*, 6(2), 1-19.

## Examples

``` r
s <- sim_kecc(n_providers = 20, type = "survival", seed = 1)
ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", time = "stop", start = "start",
              event = "event", cutoff = 0, message = FALSE)
fit <- fit_cox(ps)
head(standardize(fit))
#>   provider observed  expected     ratio person_time
#> 1     P001        6  9.113819 0.6583409    17.89913
#> 2     P002       21 24.908704 0.8430788    51.13512
#> 3     P003        8  4.686770 1.7069324    12.64318
#> 4     P004       28 15.174895 1.8451528    35.70869
#> 5     P005       22 15.468878 1.4222104    33.76398
#> 6     P006       21 24.721026 0.8494793    46.92523
head(standardize(fit, type = "direct"))
#>   provider observed expected     ratio n_pop
#> 1     P001      336 199.1986 0.5928530  2679
#> 2     P002      336 288.2112 0.8577714  2679
#> 3     P003      336 531.4106 1.5815790  2679
#> 4     P004      336 596.6812 1.7758369  2679
#> 5     P005      336 475.9018 1.4163744  2679
#> 6     P006      336 277.3695 0.8255045  2679
```
