# Confidence intervals for the provider effects of a random-effect fit

Wald intervals \\\hat\alpha_i \pm q \sqrt{v_i}\\ from the predicted
random effects and their conditional variances.

## Usage

``` r
# S3 method for class 'logistic_re'
confint(
  object,
  parm = NULL,
  level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  ...
)

# S3 method for class 'linear_re'
confint(
  object,
  parm = NULL,
  level = 0.95,
  alternative = c("two.sided", "greater", "less"),
  ...
)
```

## Arguments

- object:

  a model from
  [`fit_logistic_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md),
  [`fit_logistic_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md),
  [`fit_linear_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md)
  or
  [`fit_linear_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md).

- parm:

  optional subset of provider identifiers.

- level:

  confidence level. Default 0.95.

- alternative:

  `"two.sided"` (default), `"greater"` (lower limit only) or `"less"`
  (upper limit only).

- ...:

  unused.

## Value

A data frame with `provider`, `estimate`, `lower`, `upper`.

## Examples

``` r
# \donttest{
d <- sim_kecc(n_providers = 10, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
confint(fit_logistic_re(pd))
#>    provider    estimate       lower      upper
#> 1      P001  0.26238127 -0.18710243  0.7118650
#> 2      P002 -0.07478071 -0.38185670  0.2322953
#> 3      P003  0.44638468 -0.04169310  0.9344625
#> 4      P004 -0.47936728 -0.83778236 -0.1209522
#> 5      P005 -0.08037700 -0.45123105  0.2904771
#> 6      P006  0.18654418 -0.10164380  0.4747322
#> 7      P007  0.29734210 -0.08344507  0.6781293
#> 8      P008 -0.83039498 -1.26775277 -0.3930372
#> 9      P009  0.16725662 -0.15561030  0.4901235
#> 10     P010  0.15916380 -0.19839607  0.5167237
# }
```
