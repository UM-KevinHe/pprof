# pprof

Provider profiling in R, organised by analysis step rather than by quality
measure. A measure such as a standardized readmission or fistula ratio is
assembled from the same components: prepare the data, fit a risk-adjustment
model with provider effects, compute standardized measures, test the
provider effects against a reference distribution, and plot.

This is the development branch of version 2.0.0, which supersedes `pprof`
1.0.3, `pprofR` and `grplasso`. The package is being rebuilt component by
component; see `NEWS.md` for what exists today.

## Components

| Step | Functions |
|---|---|
| Data | `pp_data()`, `pp_check()`, `sim_kecc()` |
| Models | `fit_logistic_fe()` (SerBIN/BAN, Firth), `fit_linear_fe()`, `fit_logistic_re()`, `fit_logistic_cre()`, `fit_linear_re()`, `fit_linear_cre()` |
| Standardization | `standardize()` |
| Tests | `test_providers()`, `test_covariates()`, `confint()`, `test_measure()` |
| Reference distributions | (next) |
| Selection | `select_lasso()`, `select_group_lasso()`, `cv_lasso()`, `cv_group_lasso()` |
| Reliability, plots | (next) |

## Installation

```r
remotes::install_github("UM-KevinHe/pprof", ref = "Yubo")
```

Requires a C++17 compiler; OpenMP is optional.
