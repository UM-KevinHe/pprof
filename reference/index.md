# Package index

## Data

Prepare and check the data every component reads.

- [`pp_data()`](https://um-kevinhe.github.io/pprof/reference/pp_data.md)
  : Prepare data for a provider model
- [`pp_check()`](https://um-kevinhe.github.io/pprof/reference/pp_check.md)
  : Data checks for a provider model

## Risk adjustment

One model with one effect per provider (fixed, random or correlated
random), its coefficients, and the selection of risk factors with the
provider effects unpenalized.

- [`fit_logistic_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_logistic_fe.md)
  : Fit the logistic model with fixed provider effects
- [`fit_linear_fe()`](https://um-kevinhe.github.io/pprof/reference/fit_linear_fe.md)
  : Fit the linear model with fixed provider effects
- [`fit_logistic_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md)
  [`fit_logistic_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md)
  [`fit_linear_re()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md)
  [`fit_linear_cre()`](https://um-kevinhe.github.io/pprof/reference/fit_mixed.md)
  : Fit random-effect and correlated-random-effect provider models
- [`fit_cox()`](https://um-kevinhe.github.io/pprof/reference/fit_cox.md)
  : Cox model with provider strata or a pooled baseline
- [`vcov(`*`<logistic_fe>`*`)`](https://um-kevinhe.github.io/pprof/reference/vcov.logistic_fe.md)
  : Model-based variance of a logistic fixed-effect fit
- [`vcov(`*`<linear_fe>`*`)`](https://um-kevinhe.github.io/pprof/reference/vcov.linear_fe.md)
  : Model-based variance of a linear fixed-effect fit
- [`vcov(`*`<cox>`*`)`](https://um-kevinhe.github.io/pprof/reference/vcov.cox.md)
  : Covariance of the covariate effects of a Cox fit
- [`gamma_variance()`](https://um-kevinhe.github.io/pprof/reference/gamma_variance.md)
  : Variance of the provider effects for a Wald test

### Tests of the risk factors

- [`summary(`*`<logistic_fe>`*`)`](https://um-kevinhe.github.io/pprof/reference/summary.logistic_fe.md)
  : Tests of the covariate coefficients of a logistic fixed-effect fit
- [`summary(`*`<linear_fe>`*`)`](https://um-kevinhe.github.io/pprof/reference/summary.linear_fe.md)
  : Tests of the covariate coefficients of a linear fixed-effect fit
- [`summary(`*`<logistic_re>`*`)`](https://um-kevinhe.github.io/pprof/reference/summary_mixed.md)
  [`summary(`*`<linear_re>`*`)`](https://um-kevinhe.github.io/pprof/reference/summary_mixed.md)
  : Tests of the fixed effects of a random-effect fit
- [`summary(`*`<cox>`*`)`](https://um-kevinhe.github.io/pprof/reference/summary.cox.md)
  : Tests of the covariate coefficients of a Cox fit

### Variable selection

- [`var_select()`](https://um-kevinhe.github.io/pprof/reference/var_select.md)
  : Variable selection with unpenalized provider effects
- [`cv()`](https://um-kevinhe.github.io/pprof/reference/cv.md) :
  Cross-validation of a selection path
- [`coef(`*`<pp_lasso>`*`)`](https://um-kevinhe.github.io/pprof/reference/coef.pp_lasso.md)
  [`predict(`*`<pp_lasso>`*`)`](https://um-kevinhe.github.io/pprof/reference/coef.pp_lasso.md)
  [`selected_covariates()`](https://um-kevinhe.github.io/pprof/reference/coef.pp_lasso.md)
  : Coefficients and predictions along a selection path

## Standardized measures

- [`standardize()`](https://um-kevinhe.github.io/pprof/reference/standardize.md)
  : Standardized measures of a fitted provider model
- [`standardize(`*`<logistic_fe>`*`)`](https://um-kevinhe.github.io/pprof/reference/standardize.logistic_fe.md)
  : Standardized ratios and rates of a logistic fixed-effect fit
- [`standardize(`*`<linear_fe>`*`)`](https://um-kevinhe.github.io/pprof/reference/standardize.linear_fe.md)
  : Standardized differences of a linear fixed-effect fit
- [`standardize(`*`<logistic_re>`*`)`](https://um-kevinhe.github.io/pprof/reference/standardize_mixed.md)
  [`standardize(`*`<linear_re>`*`)`](https://um-kevinhe.github.io/pprof/reference/standardize_mixed.md)
  : Standardized measures of a random-effect fit
- [`standardize(`*`<cox>`*`)`](https://um-kevinhe.github.io/pprof/reference/standardize.cox.md)
  : Standardized ratios of a Cox fit
- [`direct_expectation()`](https://um-kevinhe.github.io/pprof/reference/direct_expectation.md)
  : Direct standardization over an arbitrary standard population

## Testing

One test per provider, intervals, the reference distribution the
statistics are judged against, and the test of a standardized measure.

- [`test()`](https://um-kevinhe.github.io/pprof/reference/test.md) :
  Tests of the provider effects
- [`test(`*`<logistic_fe>`*`)`](https://um-kevinhe.github.io/pprof/reference/test.logistic_fe.md)
  : Tests of the provider effects of a logistic fixed-effect fit
- [`test(`*`<linear_fe>`*`)`](https://um-kevinhe.github.io/pprof/reference/test.linear_fe.md)
  : Tests of the provider effects of a linear fixed-effect fit
- [`test(`*`<logistic_re>`*`)`](https://um-kevinhe.github.io/pprof/reference/test_mixed.md)
  [`test(`*`<linear_re>`*`)`](https://um-kevinhe.github.io/pprof/reference/test_mixed.md)
  : Tests of the provider effects of a random-effect fit
- [`test(`*`<cox>`*`)`](https://um-kevinhe.github.io/pprof/reference/test.cox.md)
  : Tests of the standardized ratios of a Cox fit
- [`confint(`*`<logistic_fe>`*`)`](https://um-kevinhe.github.io/pprof/reference/confint.logistic_fe.md)
  : Confidence intervals for the provider effects of a logistic
  fixed-effect fit
- [`confint(`*`<linear_fe>`*`)`](https://um-kevinhe.github.io/pprof/reference/confint.linear_fe.md)
  : Confidence intervals for the provider effects of a linear
  fixed-effect fit
- [`confint(`*`<logistic_re>`*`)`](https://um-kevinhe.github.io/pprof/reference/confint_mixed.md)
  [`confint(`*`<linear_re>`*`)`](https://um-kevinhe.github.io/pprof/reference/confint_mixed.md)
  : Confidence intervals for the provider effects of a random-effect fit
- [`confint(`*`<cox>`*`)`](https://um-kevinhe.github.io/pprof/reference/confint.cox.md)
  : Confidence limits of the standardized ratios of a Cox fit
- [`test_measure()`](https://um-kevinhe.github.io/pprof/reference/test_measure.md)
  : Test standardized measures against a reference value

### Reference distributions

Theoretical, fixed and empirical nulls, and calibration of a test with
one.

- [`null_model()`](https://um-kevinhe.github.io/pprof/reference/null_models.md)
  : Null distributions for provider statistics
- [`assign_groups()`](https://um-kevinhe.github.io/pprof/reference/assign_groups.md)
  : Group providers by size
- [`calibrate()`](https://um-kevinhe.github.io/pprof/reference/calibrate.md)
  : Calibrate a test result with a null model
- [`adjust_z()`](https://um-kevinhe.github.io/pprof/reference/adjust_z.md)
  : Adjust statistics with a null model

## Reliability

- [`iur()`](https://um-kevinhe.github.io/pprof/reference/iur.md)
  [`ratio_measure()`](https://um-kevinhe.github.io/pprof/reference/iur.md)
  : Inter-unit reliability of a provider measure
- [`iur_curve()`](https://um-kevinhe.github.io/pprof/reference/iur_curve.md)
  : Reliability at representative provider sizes
- [`iur_stratified()`](https://um-kevinhe.github.io/pprof/reference/iur_stratified.md)
  : Reliability within strata of providers

## Plots

Caterpillar, funnel, observed-expected, forest, null-calibration,
reliability, selection-path and flag-stability displays.

- [`plot(`*`<pp_test>`*`)`](https://um-kevinhe.github.io/pprof/reference/plots.md)
  [`plot(`*`<pp_measure>`*`)`](https://um-kevinhe.github.io/pprof/reference/plots.md)
  [`plot(`*`<pp_summary>`*`)`](https://um-kevinhe.github.io/pprof/reference/plots.md)
  [`plot(`*`<pp_null>`*`)`](https://um-kevinhe.github.io/pprof/reference/plots.md)
  [`plot(`*`<pp_iur>`*`)`](https://um-kevinhe.github.io/pprof/reference/plots.md)
  [`plot(`*`<pp_lasso>`*`)`](https://um-kevinhe.github.io/pprof/reference/plots.md)
  [`plot(`*`<pp_cv>`*`)`](https://um-kevinhe.github.io/pprof/reference/plots.md)
  [`plot(`*`<pp_stability>`*`)`](https://um-kevinhe.github.io/pprof/reference/plots.md)
  : Plots of profiling results
- [`flag_stability()`](https://um-kevinhe.github.io/pprof/reference/flag_stability.md)
  : Flags of several tests side by side

## Simulation

- [`sim_kecc()`](https://um-kevinhe.github.io/pprof/reference/sim_kecc.md)
  : Simulate data shaped like a KECC quality measure

## Package

- [`pprof`](https://um-kevinhe.github.io/pprof/reference/pprof-package.md)
  [`pprof-package`](https://um-kevinhe.github.io/pprof/reference/pprof-package.md)
  : pprof: provider profiling by components
