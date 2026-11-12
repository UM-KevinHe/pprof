# Standardized ratios and rates of a discrete-time fixed-effect fit

Indirect or direct standardization of every provider's events over the
person-periods at risk, with the observed and expected events and
optional confidence intervals mapped from the limits of the provider
effects.

## Usage

``` r
# S3 method for class 'discrete_fe'
standardize(
  fit,
  type = c("indirect", "direct"),
  null = "median",
  n_pop = NULL,
  providers = NULL,
  gamma_se = NULL,
  level = NULL,
  ci = c("score", "robust_score", "wald", "robust_wald"),
  alternative = c("two.sided", "greater", "less"),
  threads = 1,
  ...
)
```

## Arguments

- fit:

  a model fitted by
  [`fit_discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_discrete_fe.md).

- type:

  `"indirect"` (default) or `"direct"`.

- null:

  for indirect standardization, the norm of the provider effect:
  `"median"` (default), `"mean"` (weighted by person-periods) or a
  number.

- n_pop:

  denominator of the rates. Defaults to the number of person-periods at
  risk.

- providers:

  optional subset of provider identifiers.

- gamma_se:

  optional standard errors of the provider effects, named by provider;
  with them the direct rate's standard error is returned.

- level:

  if given, confidence limits of the ratio and the rate are added,
  mapped from the limits of \\\gamma_i\\ obtained with `ci`.

- ci:

  the interval used for the limits, see
  [`confint.discrete_fe()`](https://um-kevinhe.github.io/pprof/reference/confint.discrete_fe.md).

- alternative:

  `"two.sided"` (default), `"greater"` or `"less"`.

- threads:

  OpenMP threads for direct standardization.

- ...:

  unused.

## Value

A data frame with the columns of
[`standardize.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/standardize.logistic_fe.md).

## Details

With \\H_j(g)\\ the expected events of subject \\j\\ over its periods at
risk under the provider effect \\g\\:

- indirect: \\E_i = \sum\_{j \in i} H_j(\gamma_0)\\, ratio \\O_i /
  E_i\\;

- direct: \\E_i = \sum_j H_j(\hat\gamma_i)\\ over every subject, ratio
  \\E_i / \sum O\\.

In both cases rate = ratio \\\times \sum O / n\_{pop}\\, the events per
person-period. These are the measures of
[`standardize.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/standardize.logistic_fe.md)
on the person-period data.
