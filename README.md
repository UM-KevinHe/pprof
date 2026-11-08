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
| Models | (next) |
| Standardization | `standardize()` |
| Tests | `test_providers()`, `test_covariates()` |
| Reference distributions | (next) |
| Reliability, selection, plots | (next) |

## Installation

```r
remotes::install_github("UM-KevinHe/pprof", ref = "Yubo")
```

Requires a C++17 compiler; OpenMP is optional.
