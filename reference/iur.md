# Inter-unit reliability of a provider measure

The fraction of the variation of a provider-level measure that reflects
true differences between providers rather than within-provider sampling
noise, \\\mathrm{IUR} = \sigma^2_b / (\sigma^2_b + \sigma^2_w / n')\\,
from an ANOVA-type decomposition of the measure's variance (He et al.
2019). Three estimators, chosen by `method`:

- `"direct"` uses provider-level estimates and their standard errors
  (for instance a standardized rate and the `se_rate` of
  [`standardize()`](https://um-kevinhe.github.io/pprof/reference/standardize.md));

- `"bootstrap"` resamples the rows of every provider with replacement
  and takes the variance of the resampled measures as the
  within-provider variance;

- `"split_half"` splits every provider's rows into two random halves and
  converts the agreement of the two halves' measures across providers
  into a reliability by the Spearman-Brown formula \\2r / (1 + r)\\; six
  agreement coefficients are reported (weighted kappa, Kendall's tau and
  Spearman's rho of the measures cut into ordered categories, and
  Pearson's r, Kendall's tau and Spearman's rho of the measures
  themselves).

## Usage

``` r
iur(
  method = c("direct", "bootstrap", "split_half"),
  estimate = NULL,
  se = NULL,
  size = NULL,
  obs = NULL,
  expected = NULL,
  provider = NULL,
  n_boot = 100,
  n_iter = 10,
  category_probs = c(0.1, 0.3, 0.7, 0.9),
  measure = NULL,
  seed = NULL
)

ratio_measure(obs, expected, provider)
```

## Arguments

- method:

  `"direct"`, `"bootstrap"` or `"split_half"`.

- estimate, se:

  provider-level estimates and standard errors (direct method).

- size:

  number of observations per provider.

- obs, expected:

  observed and expected outcome per row (patient; bootstrap and
  split-half methods).

- provider:

  provider identifier per row.

- n_boot:

  number of bootstrap replicates. Default 100.

- n_iter:

  number of random splits. Default 10.

- category_probs:

  quantile probabilities that cut the measures into ordered categories
  for the categorical coefficients.

- measure:

  a function `(obs, expected, provider)` returning one value per
  provider in sorted provider order; by default the ratio of sums \\\sum
  O / \sum E\\.

- seed:

  optional seed; the session's random-number state is restored
  afterwards.

## Value

An object of class `pp_iur`: a list with `method`, `iur` (for the
split-half method a named vector of the six coefficients'
reliabilities), and for the direct and bootstrap methods `s2_between`,
`s2_within`, `n_prime`, `n_providers`, `providers`, `sizes`, `measure`,
`iur_providers` (each provider's reliability at its own size) and, for
the bootstrap, the `replicates` matrix; for the split-half method
`iterations` (one row per split) and `excluded` (providers with fewer
than two rows).

## Details

With sizes \\n_i\\, within variances \\v_i\\ of the provider-level
measure and measures \\T_i\\: the pooled within variance \\s^2\_{w} =
\sum (n_i - 1) v_i / \sum (n_i - 1)\\, the effective size \\n' = (N -
\sum n_i^2 / N) / (m - 1)\\, the total variance \\s^2_t = \sum n_i
(T_i - \bar T)^2 / (n'(m - 1))\\ with \\\bar T\\ the size-weighted mean,
the between variance \\s^2_b = s^2_t - s^2_w\\, and \\\mathrm{IUR} =
s^2_b / s^2_t\\. The reported `s2_within` is \\n' s^2_w\\, the
record-level within variance, so that the reliability of a provider of
size \\n\\ is \\s^2_b / (s^2_b + s2\\within / n)\\ (see
[`iur_curve()`](https://um-kevinhe.github.io/pprof/reference/iur_curve.md)).

## References

He K, Kalbfleisch JD, Yang Y, Li Y (2019). Inter-unit reliability for
nonlinear models. *Statistics in Medicine*, 38(5), 844-854.

## See also

[`iur_curve()`](https://um-kevinhe.github.io/pprof/reference/iur_curve.md),
[`iur_stratified()`](https://um-kevinhe.github.io/pprof/reference/iur_stratified.md)

## Examples

``` r
d <- sim_kecc(n_providers = 30, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_fe(pd)
expected <- stats::median(coef(fit)$gamma) + fit$linear_predictor
expected <- stats::plogis(expected)
boot <- iur("bootstrap", obs = fit$data$y, expected = expected,
            provider = fit$data$provider, n_boot = 50, seed = 1)
boot
#> Bootstrap inter-unit reliability: 0.8444 (30 providers)
#>   between variance 0.101, within variance 2.682, effective size 144.1
iur_curve(boot)
#>     size_group      size       iur
#> 1          min  51.00000 0.6575734
#> 2   quantile 1  61.66667 0.6989740
#> 3   quantile 2  80.66667 0.7523153
#> 4   quantile 3 108.75000 0.8037226
#> 5   quantile 4 114.00000 0.8110538
#> 6   quantile 5 128.66667 0.8289067
#> 7   quantile 6 159.33333 0.8571321
#> 8   quantile 7 168.66667 0.8639624
#> 9   quantile 8 198.66667 0.8820827
#> 10  quantile 9 208.00000 0.8867747
#> 11 quantile 10 221.00000 0.8927205
#> 12         max 231.00000 0.8968858
```
