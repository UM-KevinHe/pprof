# Null distributions for provider statistics

A null model gives, for every provider, the mean and standard deviation
of the normal distribution its z statistic is referred to: \\z\_{adj} =
(z - \mathrm{null\\mean}) / \mathrm{null\\sd}\\. Three kinds exist: the
theoretical null \\N(0, 1)\\, a fixed null with values given by the
caller (for instance a known scale inflation), and the empirical null
estimated robustly from the statistics themselves, overall or within
groups of providers of similar size. Any of them can be handed to
[`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md),
[`test()`](https://um-kevinhe.github.io/pprof/reference/test.md) or
[`test_measure()`](https://um-kevinhe.github.io/pprof/reference/test_measure.md).

## Usage

``` r
null_model(
  kind = c("empirical", "theoretical", "fixed"),
  z = NULL,
  size = NULL,
  groups = NULL,
  n_groups = 4,
  grouping = c("rank", "quantile", "cut"),
  estimator = c("huber", "bisquare", "mm"),
  tuning = NULL,
  maxit = 20,
  acc = 1e-04,
  common_mean = NULL,
  fit_mask = NULL,
  min_group_size = 3,
  small_group = c("error", "theoretical"),
  mean = 0,
  sd = 1
)
```

## Arguments

- kind:

  `"empirical"` (default), `"theoretical"` or `"fixed"`.

- z:

  for the empirical null, the provider statistics: a numeric vector,
  optionally named by provider, or a test result from
  [`test()`](https://um-kevinhe.github.io/pprof/reference/test.md) or
  [`test_measure()`](https://um-kevinhe.github.io/pprof/reference/test_measure.md),
  whose statistic, provider identifiers and (for
  [`test()`](https://um-kevinhe.github.io/pprof/reference/test.md))
  sizes are used.

- size:

  one size per provider used to form the groups (patients, say); taken
  from a test result when `z` is one. Ignored when `groups` is given.

- groups:

  explicit group label per provider; missing labels get no parameters.

- n_groups:

  number of size groups. Default 4; 1 fits one overall null.

- grouping:

  `"rank"` (default; the providers sorted by size are dealt into
  `n_groups` consecutive blocks of `ceiling(n / n_groups)`, the last
  block holding the remainder, as the KECC production code for the PPPW
  and SFR measures and pprof_py do), `"quantile"` (groups by the
  empirical quantiles of `size`, a size equal to a break going to the
  lower group, as the production SMR and SHR code and `EmpiNull` do) or
  `"cut"` (the sorted providers cut by
  [`cut()`](https://rdrr.io/r/base/cut.html) into `n_groups` intervals
  of equal length, the `Group_EN` of the KECC measures paper).

- estimator:

  robust location-scale estimator of each group, all fitted with
  [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html) on an
  intercept-only model: `"huber"` (default; Huber's psi with the MAD
  scale, the production choice), `"bisquare"` (Tukey's biweight with the
  MAD scale) or `"mm"` (MM-estimation: an S-estimate of scale followed
  by bisquare steps).

- tuning:

  tuning constant of the psi function; `MASS` defaults 1.345 (Huber) and
  4.685 (bisquare).

- maxit, acc:

  iteration limit and convergence tolerance of
  [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html) (defaults 20
  and 1e-4; `EmpiNull`'s robust calibration uses 1000 and 1e-8 with the
  bisquare).

- common_mean:

  `NULL` (default; each group's own fitted location), `TRUE` (the mean
  of all eligible statistics) or a number: replaces the groups'
  locations, the scales are kept as fitted.

- fit_mask:

  logical, one per provider: the providers that enter the fits (`TRUE`);
  the others are still calibrated with their group's parameters.

- min_group_size:

  the smallest number of eligible statistics a group may have. Default
  3.

- small_group:

  `"error"` (default) or `"theoretical"` (a group below `min_group_size`
  gets the theoretical null, with a warning).

- mean, sd:

  mean and standard deviation of a fixed null.

## Value

An object of class `pp_null`: a list with `kind` (`"theoretical"`,
`"fixed"` or `"empirical"`), `mean`, `sd` and `group` (one per provider
for the empirical null, named when the statistics were; single values
otherwise), `diagnostics` (one row per group: `group`, `n_providers`,
`n_fitted`, `min_size`, `median_size`, `max_size`, `mean_z`, `sd_z`,
`null_mean`, `null_sd`, `converged`, `n_iter`, `note`) and `settings`.

## References

Efron B (2004). Large-scale simultaneous hypothesis testing: the choice
of a null hypothesis. *Journal of the American Statistical Association*,
99(465), 96-104.

Kalbfleisch JD, Wolfe RA (2013). On monitoring outcomes of medical
providers. *Statistics in Biosciences*, 5, 286-302.

## See also

[`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md),
[`adjust_z()`](https://um-kevinhe.github.io/pprof/reference/adjust_z.md)

## Examples

``` r
d <- sim_kecc(n_providers = 40, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_fe(pd)
tst <- test(fit, statistic = "score", variance = "robust")
nul <- null_model("empirical", tst, n_groups = 2)
nul
#> Empirical null: 2 group(s) by rank, huber estimator
#>  group n_providers n_fitted min_size max_size null_mean null_sd converged
#>      1          20       20       15       39  -0.04996   2.523      TRUE
#>      2          20       20       42       60   0.11395   2.741      TRUE
head(calibrate(tst, nul))
#> Provider tests: score statistic, robust variance (modified); null median = -0.965; two.sided, level 0.95
#>   6 providers: 0 above, 0 below, 6 not different, 0 not tested
#>  provider estimate se null_value statistic z_adjusted null_mean null_sd
#>      P001  -0.7418 NA     -0.965    0.4921     0.2148  -0.04996   2.523
#>      P002  -0.3816 NA     -0.965    2.8172     0.9860   0.11395   2.741
#>      P003  -1.7131 NA     -0.965   -2.1481    -0.8315  -0.04996   2.523
#>      P004  -1.9726 NA     -0.965   -3.7949    -1.4258   0.11395   2.741
#>      P005  -1.5766 NA     -0.965   -2.8597    -1.1135  -0.04996   2.523
#>      P006  -0.6088 NA     -0.965    1.8824     0.6451   0.11395   2.741
#>  null_group      p flag  lower   upper observed expected variance
#>           1 0.8299    0   -Inf     Inf       20    17.15    33.51
#>           2 0.3241    0   -Inf  1.1226       74    51.26    65.16
#>           1 0.4057    0   -Inf     Inf       10    17.22    11.28
#>           2 0.1539    0   -Inf -0.1752       16    35.57    26.61
#>           1 0.2655    0   -Inf  0.3436       23    36.29    21.59
#>           2 0.5189    0 -4.818  0.7367       71    56.48    59.54
```
