# pprof 2.0.0 (development)

A rewrite of `pprof` 1.0.3 that absorbs `pprofR` 0.2.0 and `grplasso` 1.0.
The package is organised by analysis step, and every function takes a data
frame with column names.

## Data

* `pp_data()` replaces `fe_data_prep()` and the data handling inside the
  1.0.3 fitters: one prepared object (rows sorted by provider and patient in
  the C locale, patient rows kept in input order, provider sizes in patients
  or rows, `cutoff` with `>=`, flags for providers with no or all events).
* `pp_check()` replaces `data_check()`: no `caret` or `olsrr`, thresholds as
  arguments, metrics returned.
* `sim_kecc()` simulates data shaped like the KECC measures (patient-month
  binary or binomial days outcomes, patient random effect).

## Logistic fixed-effect model (the `pprofR` engine under the unified interface)

* `fit_logistic_fe()` replaces `logis_fe()`: binary or binomial outcomes, SerBIN or BAN.
* `test_providers()` replaces `test_fe()`: `statistic` (exact, score, Wald) x `variance` (naive,
  robust, stabilized) with one fixed set of output columns; the naive Wald statistic uses the
  case-mix variance (`gamma_variance()`), which does not depend on the covariate origin; the
  restricted score test of 1.0.3 is back as `score_type = "restricted"`.
* `test_covariates()` replaces `summary_fe_covar()`; the likelihood-ratio refit respects the fit's
  own `bound`.
* `confint()` gives intervals for the provider effects; `standardize()` replaces `SR_output()` and
  maps those limits to the ratios and rates; `direct_expectation()` replaces `direct_exp()`.
* `test_measure()` replaces `test_SM()` and adds the log scale.

## Linear, random-effect and correlated random-effect models (rewritten from 1.0.3)

* `pp_data()` detects or takes an `outcome` type (binary, binomial, continuous); `sim_kecc()`
  simulates continuous outcomes.
* `fit_linear_fe()` replaces `linear_fe()`: the within-provider regression instead of a dense
  block-diagonal centring matrix; `gamma_variance()` offers the case-mix variance
  (1.0.3's "simplified", exactly `sigma^2 / n_i`) and the "full" one; the provider and covariate
  tests use the t distribution with `n - m - p` degrees of freedom for every model-based
  statistic (1.0.3 mixed normal and t references); a patient-clustered robust variance is new.
* `fit_logistic_fe(penalty = "firth")` replaces `logis_firth()` and extends it to binomial
  trials; the object is a `logistic_fe` fit, so every test and standardization applies. The kernel
  follows ParaFLR (UM-KevinHe/hpc-flr): per-provider work on independent threads, with the p x p
  reductions as whole-data matrix products so that the result does not depend on the thread count.
* `fit_logistic_re()`, `fit_logistic_cre()`, `fit_linear_re()`, `fit_linear_cre()` replace
  `logis_re()`, `logis_cre()`, `linear_re()`, `linear_cre()` (lme4 wrappers): the random effects
  are `coefficients$alpha`, their conditional variances `post_var`; the CRE decomposition names
  its columns `<x>_within` and `<x>_between` and defaults to every covariate; the tests of the
  fixed effects use `2 * pnorm(-|z|)` (1.0.3 computed `2 * (1 - pnorm(z))`, wrong for negative
  statistics); standardization offers the observed (default) or the predicted numerator.
* The direct standardized difference of the linear models is relative to the population (the
  row-weighted mean of the effects), as the direct ratio of the logistic models always was;
  1.0.3 used the median for the linear model only.

## Variable selection (rewritten from grplasso 1.0)

* `select_lasso()` and `select_group_lasso()` replace `pp.lasso()` and `grp.lasso()`: the same
  objective and coordinate-descent paths, now with binomial trials, a deterministic update order
  (grplasso permuted the coordinates at random inside the C++ code) and without the data race
  in its OpenMP loop; `penalize` and `group` name the covariates instead of positional vectors.
* `cv_lasso()` and `cv_group_lasso()` replace `cv.pp.lasso()` and `cv.grp.lasso()`: folds are
  formed within each provider over patients (grplasso stratified rows on the outcome and could
  leave a provider out of a training set, which crashed its predict method), the deviance is
  averaged by fold with a fold-based standard error, and both `lambda_min` and `lambda_1se` are
  reported.
* `coef()`, `predict()` and `selected_covariates()` work along the path; `predict()` takes a
  data frame and handles unseen providers explicitly (grplasso indexed the effects by position
  for numeric identifiers).
* Discrete-time survival and stratified Cox selection move to the survival components.

## Still to come (in this order)

* Empirical null, inter-unit reliability, Cox O/E components.
* Plots.
