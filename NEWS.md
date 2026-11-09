# pprof 2.0.6

* Score tests carry confidence limits: `test()` inverts the modified score
  statistic (naive, robust or stabilized variance) for every provider, so the
  interval excludes the null value exactly when the provider is flagged, and
  the caterpillar plot of a score test shows intervals. `calibrate()` inverts
  the calibrated statistic the same way. `confint()` gains
  `type = "robust_score"` and `"stabilized_score"`. `interval = FALSE` in
  `test()` skips the root finding.

# pprof 2.0.5

* Get started: the analysis tests with the robust Wald statistic, so the
  caterpillar plot carries intervals; the plots section shows one plain test
  (caterpillar, funnel, observed against expected) and the empirical-null
  figure moves to the Reference distributions article.

# pprof 2.0.4

* README: the ten-line example and its figures use the plain score test under the
  theoretical null.

# pprof 2.0.3

* The funnel plot of a score test (`plot(test, type = "funnel")`) is drawn:
  observed over expected events against the precision $E^2 / V$, with $V$
  the variance the statistic was referred to (naive, robust or stabilized),
  and the limits $1 + (m \pm c\, s) / \sqrt{E^2 / V}$ of its null (one set
  per group after `calibrate()`), so that a provider lies outside the limits
  exactly when it is flagged. `test()` results of the score and exact
  statistics carry `observed`, `expected` and `variance` for it.
* `plot(measure, type = "funnel", test = )` draws that test's funnel instead
  of normal limits with the test's colours, which could show many providers
  outside the limits and few flagged; a Wald-type test is refused there.
* One-sided Wald funnels draw the limit on the flagged side.
* README funnel figure redrawn with the test's own limits.

# pprof 2.0.2

* README and site home page shortened to the essentials, with the plots of the
  Get started analysis (`man/figures/`).

# pprof 2.0.1

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
* `test()` replaces `test_fe()`: `statistic` (exact, score, Wald) x `variance` (naive,
  robust, stabilized) with one fixed set of output columns; the naive Wald statistic uses the
  case-mix variance (`gamma_variance()`), which does not depend on the covariate origin; the
  restricted score test of 1.0.3 is back as `score_type = "restricted"`.
* `summary()` replaces `summary_fe_covar()`; the likelihood-ratio refit respects the fit's
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
  `test(statistic = "lr")` adds ParaFLR's likelihood-ratio test of one provider effect
  with the other parameters fixed (penalized likelihood for a Firth fit).
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

* `var_select()` (lasso or group lasso) replace `pp.lasso()` and `grp.lasso()`: the same
  objective and coordinate-descent paths, now with binomial trials, a deterministic update order
  (grplasso permuted the coordinates at random inside the C++ code) and without the data race
  in its OpenMP loop; `penalize` and `group` name the covariates instead of positional vectors.
* `cv()` replace `cv.pp.lasso()` and `cv.grp.lasso()`: folds are
  formed within each provider over patients (grplasso stratified rows on the outcome and could
  leave a provider out of a training set, which crashed its predict method), the deviance is
  averaged by fold with a fold-based standard error, and both `lambda_min` and `lambda_1se` are
  reported.
* `coef()`, `predict()` and `selected_covariates()` work along the path; `predict()` takes a
  data frame and handles unseen providers explicitly (grplasso indexed the effects by position
  for numeric identifiers).
* Discrete-time survival and stratified Cox selection move to the survival components.

## Reference distributions and reliability (new; empirical null from the KECC production code and pprof_py)

* Null models `null_model()` (theoretical, fixed or empirical): the empirical null is fitted
  with `MASS::rlm` within size groups. The default (`grouping = "rank"`: the providers sorted by size
  dealt into blocks of `ceiling(n / 4)`; `estimator = "huber"` with MASS defaults) is what the
  production PPPW and SFR code does; `grouping = "quantile"` with `estimator = "bisquare"`,
  `maxit = 1000`, `acc = 1e-8` is the production SMR code (and `EmpiNull`), `"quantile"` with the
  Huber default the production SHR code, and `grouping = "cut"` the `Group_EN` of the KECC measures
  paper. MM-estimation, tuning constants, a common mean, a fit mask and small-group handling are the
  options of pprof_py's `EmpiricalNull`.
* `calibrate()` refers a `test()` or `test_measure()` result to a null model (adjusted
  statistics, p-values, flags and inverted limits); every `test()` method and
  `test_measure()` take a `null_model`. `critical` sets the critical value of the flags and
  limits; `ci_location` centres the limits elsewhere than the null mean (`0` gives limits of
  the form estimate ± c · null_sd · se).
* Inter-unit reliability translated from pprof_py: `iur()` with the direct, bootstrap and
  split-half methods, with `iur_curve()` and `iur_stratified()`.

## Still to come (in this order)

* Cox O/E components, discrete-time survival and stratified Cox selection.
* Plots.
* `var_select()` (lasso or group lasso) on a survival outcome fit `grplasso`'s
  `Strat.cox()` (the Cox model stratified by provider, now with counting-process rows, case
  weights and offsets) and, with `model = "discrete"`, `pp.DiscSurv()` (the discrete-time hazard
  model with period and provider effects, the provider effects identified by a zero median
  instead of a reference provider), both with a cyclic update order and without the hard-coded
  column names; `cv()` cross-validate them on the partial-likelihood
  deviance per event or the person-period deviance.

## Interface

* One entry point per analysis step, with S3 methods for every model class: `test()` (provider
  effects; was `test_providers()`), `summary()` (covariate coefficients; was `test_covariates()`),
  `standardize()`, `confint()`, `calibrate()`, `coef()`, `vcov()`, `predict()`; `null_model(kind)`
  for the theoretical, fixed and empirical nulls; `iur(method)` for the three reliability
  estimators; `var_select(penalty)` for the lasso and group-lasso paths with `cv()` on the path.
  `test()` and `test_measure()` also accept a null kind by name (`null_model = "empirical"`).

## Plots (ggplot2)

* `plot()` methods, after pprof_py's presentation layer and the 1.0.3 plots: a `test()` result
  gives a caterpillar plot, a funnel plot (limits from the same test, on the precision or the
  expected-count axis, one curve per empirical-null group when calibrated) or the flag shares by
  provider-size group; a `standardize()` result the observed-against-expected scatter or an
  observed-over-expected funnel; a `summary()` result a forest plot of the coefficients; an
  empirical `null_model()` the null-calibration panels; an `iur()` result the reliability curve;
  a `var_select()` path and its `cv()` the coefficient paths and the cross-validation curve.
  `flag_stability()` tabulates the flags of several tests (scenarios) and plots them.
