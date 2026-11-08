# Fit the logistic model with fixed provider effects

Fits \\\mathrm{logit}\\ p\_{ij} = \gamma_i + Z\_{ij}^\top\beta\\ with a
binary or binomial outcome (\\Y\_{ij} \sim \mathrm{Bin}(N\_{ij},
p\_{ij})\\) by the serial blockwise inversion Newton (SerBIN) or block
ascent Newton (BAN) algorithm, on the providers that
[`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
marked as included, optionally with Firth's penalty.

## Usage

``` r
fit_logistic_fe(
  data,
  algorithm = c("SerBIN", "BAN"),
  penalty = c("none", "firth"),
  max_iter = 10000,
  tol = 1e-05,
  bound = 10,
  backtrack = TRUE,
  threads = 1,
  stop = c("beta", "relch", "ratch", "all", "or"),
  auc = FALSE,
  message = FALSE
)
```

## Arguments

- data:

  an object from
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  with a binary or binomial outcome.

- algorithm:

  `"SerBIN"` (default; Wu et al. 2022) or `"BAN"` (He et al. 2013).

- penalty:

  `"none"` (default; maximum likelihood) or `"firth"` (penalized
  likelihood with Jeffreys' prior, Firth 1993), which keeps the effects
  of providers with no or all events finite.

- max_iter:

  maximum number of iterations. Default 10,000.

- tol:

  convergence tolerance. Default 1e-5.

- bound:

  provider effects are kept within `median(gamma) +/- bound` during the
  iterations. Default 10.

- backtrack:

  use a backtracking line search (not used with the Firth penalty).
  Default `TRUE`.

- threads:

  number of OpenMP threads for the information matrix (SerBIN and
  Firth). Results do not depend on it.

- stop:

  stopping rule: `"beta"` (default; infinity norm of the change in
  beta), `"relch"` (relative change in the log-likelihood), `"ratch"`
  (change relative to the total change since the start), `"all"` (all
  three below `tol`) or `"or"` (any of the three below `tol`). With the
  Firth penalty the rules read the penalized log-likelihood.

- auc:

  also compute the area under the ROC curve (binary outcomes).

- message:

  print the iteration history.

## Value

An object of class `logistic_fe`:

- coefficients:

  a list with `gamma` (provider effects, named by provider) and `beta`
  (covariate coefficients, named by covariate).

- fitted:

  fitted probabilities of the included rows.

- linear_predictor:

  \\Z\hat\beta\\ of the included rows.

- loglik, aic, bic:

  binomial log-likelihood (without the constant \\\log\binom{N}{Y}\\)
  and the information criteria.

- penalized_loglik:

  the penalized log-likelihood (Firth only).

- auc:

  the AUC, when requested (`NA` for binomial outcomes).

- data:

  the included rows, as sorted by
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md).

- vars, sizes, cutoff:

  from
  [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  (`sizes` restricted to the included providers).

- rows:

  rows per provider, in the order of `coefficients$gamma`.

- control:

  the settings and the number of iterations.

- converged:

  whether the stopping rule was met before `max_iter`.

## Details

SerBIN updates provider effects and covariate coefficients jointly using
the block structure of the information matrix; BAN alternates between
the two blocks. Providers with no or all events are fitted (their
effects sit at the bound); the tests treat them as the data allow.

With `penalty = "firth"` the penalized log-likelihood \\\ell(\gamma,
\beta) + \frac12 \log\det I(\gamma, \beta)\\ is maximized by a joint
Newton step through the same block structure; the modified score adds
\\h\_{ij}(\frac12 - p\_{ij})\\ to each row's residual, with \\h\_{ij}\\
the row's leverage \\N\_{ij} p\_{ij}(1 - p\_{ij})\\ x\_{ij}^\top I^{-1}
x\_{ij}\\. The estimates stay finite when a provider has no or all
events. The computation follows ParaFLR (Shao and He, the `paraflr`
package in UM-KevinHe/hpc-flr): the provider is the unit of parallel
work, so the per-provider information blocks, leverages and scores are
spread over `threads`; the results do not depend on the number of
threads. `loglik`, `aic` and `bic` report the unpenalized log-likelihood
at the penalized estimate; the penalized value is `penalized_loglik`.

## References

Firth D (1993). Bias reduction of maximum likelihood estimates.
*Biometrika*, 80(1), 27-38.

Shao Y, He K. ParaFLR: parallel Firth-corrected logistic regression with
high-dimensional provider effects.
<https://github.com/UM-KevinHe/hpc-flr>

He K, Kalbfleisch JD, Li Y, Li Y (2013). Evaluating hospital readmission
rates in dialysis facilities; adjusting for hospital effects. *Lifetime
Data Analysis*, 19, 490-512.

Wu W, Yang Y, Kang J, He K (2022). Improving large-scale estimation and
inference for profiling health care providers. *Statistics in Medicine*,
41(15), 2840-2853.

## See also

[`test()`](https://um-kevinhe.github.io/pprof/reference/test.md),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`standardize()`](https://um-kevinhe.github.io/pprof/reference/standardize.md),
[`confint.logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/confint.logistic_fe.md)

## Examples

``` r
d <- sim_kecc(n_providers = 20, seed = 1)
pd <- pp_data(d, y = "y", x = c("x1", "x2", "x3"), provider = "provider",
              patient = "patient", message = FALSE)
fit <- fit_logistic_fe(pd)
fit
#> Logistic model with fixed provider effects (SerBIN)
#>   2679 rows, 20 providers, 3 covariates; binary outcome
#>   log-likelihood -1336.52; iterations 5
#>   provider effects: median -1.08, range -1.823058 to  0.006187
coef(fit)$beta
#>         x1         x2         x3 
#>  0.5942121 -0.4515167  0.3261670 
fit_firth <- fit_logistic_fe(pd, penalty = "firth")
```
