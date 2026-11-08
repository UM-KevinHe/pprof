# Fit random-effect and correlated-random-effect provider models

Wrappers around
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) (binary and
binomial outcomes) and
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) (continuous
outcomes) that fit a random provider intercept, \$\$g(E\[Y\_{ij}\]) =
\mu + \alpha_i + Z\_{ij}^\top\beta, \qquad \alpha_i \sim N(0,
\sigma\_\alpha^2),\$\$ on the providers that
[`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
marked as included. The correlated random-effect (CRE) models add the
Mundlak within-between decomposition of selected covariates: for each
decomposed covariate the provider mean (`<name>_between`) and the
deviation from it (`<name>_within`) enter the model instead of the
covariate itself, which lets the random intercept be correlated with the
covariates' provider means.

## Usage

``` r
fit_logistic_re(data, ...)

fit_logistic_cre(data, decompose = NULL, ...)

fit_linear_re(data, ...)

fit_linear_cre(data, decompose = NULL, ...)
```

## Arguments

- data:

  an object from
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md):
  a binary or binomial outcome for the logistic models, a continuous
  outcome for the linear models.

- ...:

  passed to [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html)
  or [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) (for
  instance `nAGQ` or `control`).

- decompose:

  names of the covariates to decompose (CRE models); default all
  covariates.

## Value

An object of class `logistic_re`, `logistic_cre` (which inherits from
`logistic_re`), `linear_re` or `linear_cre` (inherits from `linear_re`):

- coefficients:

  a list with `alpha` (the predicted random effects, named by provider)
  and `beta` (the fixed effects, starting with `(Intercept)`).

- variance:

  a list with `alpha` (\\\hat\sigma\_\alpha^2\\) and `beta` (the
  covariance matrix of the fixed effects).

- post_var:

  the conditional variance of each \\\hat\alpha_i\\ given the data
  (lme4's `postVar`), named by provider.

- sigma:

  the residual standard deviation (linear models).

- fitted:

  fitted means (probabilities for the logistic models) of the included
  rows, including the random effects.

- linear_predictor:

  the fixed part \\\hat\mu + Z\hat\beta\\ of the included rows.

- residuals:

  response residuals \\Y - N\hat p\\ or \\Y - \hat \mu\\.

- loglik, aic, bic:

  from the lme4 fit.

- model:

  the lme4 fit.

- data:

  the included rows, with the decomposed columns added (CRE).

- vars:

  as in
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md);
  for the CRE models `x` holds the names of the model terms and
  `decompose` the covariates decomposed.

- sizes, cutoff, rows:

  as in
  [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md).

## References

Bates D, Maechler M, Bolker B, Walker S (2015). Fitting linear
mixed-effects models using lme4. *Journal of Statistical Software*,
67(1), 1-48.

Mundlak Y (1978). On the pooling of time series and cross section data.
*Econometrica*, 46(1), 69-85.

## See also

[`test.logistic_re()`](https://um-kevinhe.github.io/pprof/reference/test_mixed.md),
[`standardize.logistic_re()`](https://um-kevinhe.github.io/pprof/reference/standardize_mixed.md),
[`confint.logistic_re()`](https://um-kevinhe.github.io/pprof/reference/confint_mixed.md),
[`summary.logistic_re()`](https://um-kevinhe.github.io/pprof/reference/summary_mixed.md)

## Examples

``` r
# \donttest{
d <- sim_kecc(n_providers = 20, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_re(pd)
fit
#> Logistic model with random provider effects
#>   2679 rows, 20 providers, 3 fixed-effect terms
#>   provider variance 0.2581
#>   log-likelihood -1364.92
fit_cre <- fit_logistic_cre(pd, decompose = "x1")
# }
```
