# Confidence limits of the standardized ratios of a Cox fit

The limits of
[`test.cox()`](https://um-kevinhe.github.io/pprof/reference/test.cox.md)
for the ratio \\O_j / E_j\\.

## Usage

``` r
# S3 method for class 'cox'
confint(
  object,
  parm = NULL,
  level = 0.95,
  statistic = c("midp", "exact", "wald"),
  dispersion = c("none", "deviance"),
  alternative = c("two.sided", "greater", "less"),
  min_expected = 0,
  null_model = NULL,
  ...
)
```

## Arguments

- object:

  a model fitted by
  [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md).

- parm:

  optional subset of provider identifiers.

- level:

  confidence level. Default 0.95.

- statistic, dispersion, alternative, min_expected, null_model:

  as in
  [`test.cox()`](https://um-kevinhe.github.io/pprof/reference/test.cox.md).

- ...:

  unused.

## Value

A data frame with `provider`, `estimate` (the ratio), `lower`, `upper`.

## Examples

``` r
s <- sim_kecc(n_providers = 20, type = "survival", seed = 1)
ps <- pp_data(s, x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", time = "stop", start = "start",
              event = "event", cutoff = 0, message = FALSE)
confint(fit_cox(ps))
#>    provider  estimate     lower     upper
#> 1      P001 0.6583409 0.2668281 1.3692670
#> 2      P002 0.8430788 0.5358011 1.2667455
#> 3      P003 1.7069324 0.7927412 3.2413570
#> 4      P004 1.8451528 1.2502304 2.6310420
#> 5      P005 1.4222104 0.9139247 2.1178970
#> 6      P006 0.8494793 0.5398688 1.2763624
#> 7      P007 0.9018641 0.4742492 1.5675270
#> 8      P008 1.1453681 0.6780111 1.8203659
#> 9      P009 2.5061235 1.7728057 3.4468898
#> 10     P010 0.4958119 0.2607248 0.8617690
#> 11     P011 1.2076866 0.7270230 1.8943837
#> 12     P012 0.5667832 0.3513649 0.8686987
#> 13     P013 0.3803804 0.2115611 0.6341139
#> 14     P014 1.4846233 0.8043992 2.5239651
#> 15     P015 0.6789168 0.3310902 1.2459198
#> 16     P016 0.9053589 0.4415203 1.6614771
#> 17     P017 1.8703490 1.1071701 2.9725985
#> 18     P018 1.0294686 0.6197364 1.6148300
#> 19     P019 1.8409337 1.2562971 2.6095226
#> 20     P020 0.6217985 0.3680793 0.9882420
```
