# Test standardized measures against a reference value

Measure-level z test. With a transformation \\g\\ (identity, logit or
log) the statistic is \$\$z_i = \frac{g(M_i) - g(M\_{ref})}{g'(M_i)\\
SE_i}.\$\$ Under the null, \\z_i\\ follows \\N(\mathrm{location},
\mathrm{scale}^2)\\: `location = 0, scale = 1` is the theoretical null,
a fixed inflation such as `scale = 1.81` is a known-scale null, and the
location and scale of an empirical null can be passed per provider.

## Usage

``` r
test_measure(
  estimate,
  se,
  reference,
  transform = c("identity", "logit", "log"),
  location = 0,
  scale = 1,
  level = 0.95,
  critical = NULL,
  ci_location = location,
  alternative = c("two.sided", "greater", "less"),
  providers = NULL,
  null_model = NULL
)
```

## Arguments

- estimate:

  the standardized measures \\M_i\\, one per provider.

- se:

  their standard errors.

- reference:

  the reference value \\M\_{ref}\\ (one number or one value per
  provider), for instance the population rate.

- transform:

  `"identity"` (default), `"logit"` or `"log"`.

- location, scale:

  location and scale of the null distribution of \\z\\ (one number or
  one value per provider).

- level:

  confidence level. Default 0.95.

- critical:

  critical value; defaults to the normal quantile for `level`.

- ci_location:

  location used to centre the confidence interval; defaults to
  `location`.

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

- providers:

  optional provider identifiers for the output.

- null_model:

  optional null model (a `pp_null` object or a function of the result
  returning one), applied through
  [`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md);
  it replaces `location` and `scale`.

## Value

A data frame of class `pp_measure_test` with one row per provider:
`provider` (if given), `estimate`, `se`, `reference`, `estimate_t`,
`se_t`, `reference_t` (on the transformed scale), `z`, `p`, `flag` (`1`
above the reference, `-1` below, `0` not different, `NA` when `z` is
missing), `lower`, `upper` (measure scale) and `lower_t`, `upper_t`.

## Details

The interval inverts the calibrated test: \\g(M_i) -
(\mathrm{ci\\location} \pm \mathrm{critical} \cdot \mathrm{scale})\\
SE_i^{(g)}\\, mapped back with \\g^{-1}\\. One-sided alternatives give
one-sided intervals.

## Examples

``` r
test_measure(estimate = c(0.12, 0.30, 0.22), se = c(0.03, 0.04, 0.05),
             reference = 0.2, transform = "logit")
#>   estimate   se reference estimate_t      se_t reference_t          z
#> 1     0.12 0.03       0.2 -1.9924302 0.2840909   -1.386294 -2.1335980
#> 2     0.30 0.04       0.2 -0.8472979 0.1904762   -1.386294  2.8297316
#> 3     0.22 0.05       0.2 -1.2656664 0.2913753   -1.386294  0.4139953
#>             p flag      lower     upper   lower_t    upper_t
#> 1 0.032875700   -1 0.07247769 0.1922242 -2.549238 -1.4356222
#> 2 0.004658706    1 0.22782660 0.3836767 -1.220624 -0.4739714
#> 3 0.678877583    0 0.13743595 0.3330147 -1.836751 -0.6945813
```
