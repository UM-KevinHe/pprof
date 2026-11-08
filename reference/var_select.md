# Variable selection with unpenalized provider effects

Lasso and group-lasso paths for the logistic model with fixed provider
effects, \\\mathrm{logit}\\ p\_{ij} = \gamma_i + Z\_{ij}^\top\beta\\,
where the provider effects (and any covariate named in `penalize`'s
complement) are never penalized. The objective at penalty \\\lambda\\ is
\$\$-\frac{1}{n}\\\ell(\gamma, \beta) + \lambda \sum_j m_j \|\beta_j\|
\qquad \text{or} \qquad -\frac{1}{n}\\\ell(\gamma, \beta) + \lambda
\sum_g m_g \\\beta_g\\\_2,\$\$ solved by coordinate descent on
standardized covariates (each penalized group orthonormalized) with one
Newton step on the provider effects per cycle, warm starts along a
decreasing sequence of \\\lambda\\, and an active-set strategy.

## Usage

``` r
var_select(
  data,
  penalty = c("lasso", "group_lasso"),
  group = NULL,
  penalize = NULL,
  multiplier = NULL,
  lambda = NULL,
  n_lambda = 100,
  lambda_min_ratio = 1e-04,
  standardize = TRUE,
  bound = 10,
  tol = 1e-04,
  max_iter = 10000,
  active_set = TRUE,
  mm = FALSE,
  backtrack = FALSE,
  threads = 1,
  model = NULL
)
```

## Arguments

- data:

  an object from
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  with a binary, binomial or survival outcome.

- penalty:

  `"lasso"` (default) or `"group_lasso"`.

- group:

  for the group lasso, the group of each covariate: a vector of labels
  in the order of the covariates, or named by covariate. The label `0`,
  `NA` or `""` marks an unpenalized covariate.

- penalize:

  names of the covariates that carry a penalty; default all. The others
  enter the model unpenalized.

- multiplier:

  penalty multipliers \\m_j\\ (named by covariate) or \\m_g\\ (named by
  group); default 1 per covariate and the square root of the group size
  per group.

- lambda:

  the penalty sequence; by default `n_lambda` values from the smallest
  \\\lambda\\ with every penalized coefficient zero down to
  `lambda_min_ratio` times it, equally spaced on the log scale.

- n_lambda, lambda_min_ratio:

  length and range of the default sequence.

- standardize:

  centre and scale the covariates before penalizing (the estimates are
  returned on the original scale). Default `TRUE`.

- bound:

  the provider effects are kept within `median +/- bound`.

- tol:

  convergence tolerance on the largest change of a coefficient within a
  cycle. Default 1e-4.

- max_iter:

  maximum number of cycles per value of `lambda`.

- active_set:

  start each `lambda` from the coefficients that are nonzero and add the
  covariates that violate the optimality conditions.

- mm:

  for the lasso, use the majorization with the constant working weight
  \\N/4\\ instead of the IRLS weights (the group lasso always does).
  Default `FALSE`.

- backtrack:

  line search in the provider-effect step. Default `FALSE`.

- threads:

  OpenMP threads for the provider-effect step.

- model:

  for a survival outcome: `"cox"` (default; the Cox model stratified by
  provider, whose baselines absorb the provider effects) or `"discrete"`
  (the discrete-time hazard model, lasso only). `NULL` for the other
  outcomes.

## Value

An object of class `pp_lasso` (and `pp_group_lasso` for the group
lasso):

- beta:

  a matrix of covariate coefficients, one column per `lambda`.

- gamma:

  a matrix of provider effects, one column per `lambda` (`NULL` for the
  stratified Cox model).

- alpha, time_grid:

  the discrete-time model's period effects (one column per `lambda`) and
  the periods' times.

- model:

  `"logistic"`, `"cox"` or `"discrete"`.

- lambda, deviance, iter:

  per value of `lambda` (the deviance of the Cox model is minus twice
  its partial log-likelihood).

- n_selected:

  number of penalized covariates with a nonzero coefficient per
  `lambda`.

- penalized:

  logical, per covariate.

- group:

  the group of each covariate (group lasso).

- multiplier:

  the multipliers used.

- data, vars, sizes, rows:

  as in
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md).

## Details

Providers with fewer than the cutoff patients are excluded as in every
fitter. The path starts from the fit with the unpenalized covariates
only (or the provider-only model), whose residuals give the largest
\\\lambda\\. The deviance reported per \\\lambda\\ is \\-2\ell\\ without
the constant term.

With a survival outcome
([`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
with `time` and `event`) the paths are those of `grplasso`'s
`Strat.cox()` and `pp.DiscSurv()`:

- `model = "cox"`: the weighted partial likelihood of the Cox model with
  one baseline per provider, counting-process rows `(start, stop]`, case
  weights and offsets as in
  [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md);
  the provider effects are not parameters (`gamma` is `NULL`) and the
  working weight is 1 (Breheny and Huang's majorization for the Cox
  model). The covariates of a group are orthonormalized as for the
  logistic model.

- `model = "discrete"`: \\\mathrm{logit}\\ h\_{jk} = \alpha_k +
  \gamma_i + Z_j^\top\beta\\ for subject \\j\\ of provider \\i\\ in
  period \\k\\, the periods being the distinct values of `time`; one row
  per subject (no `start`, no weights), the subject at risk in every
  period up to its own. The period effects \\\alpha\\ and the provider
  effects \\\gamma\\ are unpenalized, \\\gamma\\ identified by a zero
  median (`grplasso` fixed the first provider at zero) and kept within
  `+/- bound`; the likelihood is computed without expanding the data.
  The IRLS weights \\\sum_k h\_{jk}(1 - h\_{jk})\\ or, with `mm = TRUE`,
  \\t_j / 4\\.

## References

Breheny P, Huang J (2015). Group descent algorithms for nonconvex
penalized linear and logistic regression models with grouped predictors.
*Statistics and Computing*, 25, 173-187.

Yang Y, Wu W, He K (2025+). Variable selection for provider profiling
with high-dimensional fixed effects (the `grplasso` package).

## See also

[`cv()`](https://um-kevinhe.github.io/pprof/reference/cv.md),
[`coef.pp_lasso()`](https://um-kevinhe.github.io/pprof/reference/coef.pp_lasso.md),
[`predict.pp_lasso()`](https://um-kevinhe.github.io/pprof/reference/coef.pp_lasso.md)

## Examples

``` r
d <- sim_kecc(n_providers = 20, n_covariates = 6, seed = 1)
pd <- pp_data(d, y = "y", x = paste0("x", 1:6), provider = "provider",
              patient = "patient", message = FALSE)
path <- var_select(pd, n_lambda = 20)
path
#> Lasso path of the logistic model with unpenalized provider effects
#>   2679 rows, 20 providers, 6 covariates (6 penalized)
#>   20 values of lambda from 0.09666 to 9.666e-06; selected covariates from 0 to 6
coef(path, lambda = path$lambda[10])$beta
#>         x1         x2         x3         x4         x5         x6 
#>  0.5670012 -0.3142008  0.2640520  0.5155683 -0.3474599  0.2350394 
grp <- var_select(pd, penalty = "group_lasso", group = c(1, 1, 2, 2, 3, 0), n_lambda = 20)
```
