# Confidence intervals for the provider effects of a linear fixed-effect fit

Wald intervals for \\\gamma_i\\ with the model-based standard error (t
reference with \\n - m - p\\ degrees of freedom) or the
patient-clustered standard error (normal reference).

## Usage

``` r
# S3 method for class 'linear_fe'
confint(
  object,
  parm = NULL,
  level = 0.95,
  type = c("wald", "robust_wald"),
  wald_variance = c("case_mix", "full"),
  alternative = c("two.sided", "greater", "less"),
  ...
)
```

## Arguments

- object:

  a model fitted by
  [`fit_linear_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_linear_fe.md).

- parm:

  optional subset of provider identifiers.

- level:

  confidence level. Default 0.95.

- type:

  `"wald"` (default) or `"robust_wald"`.

- wald_variance:

  for `type = "wald"`: `"case_mix"` (default) or `"full"`, see
  [`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md).

- alternative:

  `"two.sided"` (default), `"greater"` (lower limit only) or `"less"`
  (upper limit only).

- ...:

  unused.

## Value

A data frame with `provider`, `estimate`, `lower`, `upper`.

## Examples

``` r
d <- sim_kecc(n_providers = 10, type = "continuous", seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
confint(fit_linear_fe(pd))
#>    provider     estimate       lower       upper
#> 1      P001  0.254765408 -0.03271138  0.54224220
#> 2      P002  0.166225574  0.00192733  0.33052382
#> 3      P003  0.009955542 -0.32457736  0.34448845
#> 4      P004  0.109453499 -0.07098013  0.28988713
#> 5      P005 -0.129846923 -0.34142411  0.08173026
#> 6      P006  0.472555209  0.31107257  0.63403785
#> 7      P007  0.579205785  0.34141356  0.81699801
#> 8      P008 -1.205787618 -1.42191537 -0.98965986
#> 9      P009  0.285299393  0.09478464  0.47581414
#> 10     P010  0.007548245 -0.20143280  0.21652929
```
