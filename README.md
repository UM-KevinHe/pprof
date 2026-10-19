# pprofR

Shared modelling, testing and standardized-measure engine for the KECC quality
measures (AOH, PPPW, SFR, SRR). It replaces the copies of `Fixed_effect.cpp`,
`fe_data_prep.R`, `logis_fe.R`, `test_fe.R`, `summary_fe_covar.R`,
`SR_output.R`, `coef.logis_fe.R`, `test_aoh.R`, `summary_aoh_beta.R` and
`AOH_piehat.cpp` that each measure folder used to carry. Empirical-null
calibration stays in the measure code.

## Install

Dependencies: `Rcpp`, `RcppArmadillo`, `poibin` (R >= 4.0). A C++17 compiler is
required. OpenMP is used when available. No `~/.R/Makevars` edits are needed:
the package sets its own C++ standard and OpenMP flags.

```r
# from the source tarball (dependencies must already be installed)
install.packages(c("Rcpp", "RcppArmadillo", "poibin"))
install.packages("pprofR_0.2.0.tar.gz", repos = NULL, type = "source")

# or, resolving dependencies automatically
remotes::install_local("pprofR_0.2.0.tar.gz")
```

## Functions

| Step | Function | Notes |
|---|---|---|
| Data | `fe_data_prep(N, Y, Z, ProvID, PatientID, cutoff, check, message)` | `N` = number of trials (binomial); omit for binary |
| Model | `logis_fe(data.prep, algorithm, tol, ...)` | SerBIN (default) or BAN; stopping rules `"beta"` (default), `"relch"`, `"ratch"`, `"all"`, `"or"` |
| | `coef()`, `vcov()` | S3 methods for `logis_fe` fits |
| Provider effects (gamma) | `test_fe(fit, test = ...)` | `"exact.poisbinom"`, `"exact.bootstrap"`, `"wald"`, `"robust_wald"`, `"stabilized_wald"`, `"score"`, `"robust_score"`, `"stabilized_score"`; `alternative =` |
| Covariates (beta) | `summary_fe_covar(fit, test = ...)` | `"wald"`, `"robust_wald"`, `"lr"`, `"score"` |
| Intervals | `confint(fit, option = "gamma" / "SM", test = ...)` | provider effects, or direct / indirect ratios and rates |
| Standardized measures | `SR_output(fit, stdz, measure, ...)` | direct / indirect ratio and rate, with the direct-rate SE numerator |
| | `direct_exp(gamma, Z_beta, N)` | direct standardization over any standard population |
| Measure-level test | `test_SM(measure, se, ref, transform, location, scale)` | z-test / p / flag / CI; identity or logit scale; null N(location, scale^2) |
| Checks | `data_check(Y, Z, ProvID)` | missingness, near-zero variance, correlation, VIF |

Binary outcomes are the special case `N = 1`; every formula reduces exactly to
the binary code the measures used before.

## How the measures use it

| Measure | Model | gamma test | beta test | Standardized measure | Measure test |
|---|---|---|---|---|---|
| SFR | `logis_fe(tol = 1e-5)` | `test_fe("robust_wald")` | – | `SR_output(direct)` | `test_SM(identity)`, then EN |
| PPPW / aPPPW | `logis_fe(tol = 1e-15)` | `test_fe("robust_wald")` | `summary_fe_covar("robust_wald")` | `SR_output(direct)` (+/-17 providers added) | `test_SM(logit, scale = 1.81 or EN)` |
| AOH | `logis_fe(N = days at risk, tol = 1e-10)` | `test_fe("robust_wald")` | `summary_fe_covar("robust_wald")` | `direct_exp(N = day weights)` | `test_SM(identity, ref = pie_l)` |
| SRR | stage 1 `logis_fe(tol = 1e-5)` | (facility FE + hospital RE stages stay in SRR code) | `summary_fe_covar("wald")` | (SRR code) | (SRR code) |

The convergence tolerance differs by measure. Each measure script passes its
own value explicitly (PPPW's former copy of `logis_fe` defaulted to 1e-15).

## Behaviour notes

- Results do not depend on the number of cores. The former sourced C++ let
  Armadillo parallelize some internal sums, so its estimates changed slightly
  with the core count (up to about 1e-10 relative in the AOH validation, where
  `tol = 1e-10`).
- Validated against the former code on simulated data (`validation/run_all.R`):
  every pprofR function agrees within 1e-14, and flags and categories are
  identical in every measure pipeline. The pipelines agree within 5e-12,
  except for values that the PPPW measure code itself derives from its +17
  placeholder providers (a logit of a rate 2e-7 below 1), where the former
  code's rounding is magnified; see the validation page of the docs.
- The robust and stabilized tests need `PatientID` in `fe_data_prep()`.
  Variances cluster on (provider, patient). The stabilized tests pool one
  within-patient correlation over all patients by row position, so a
  patient's rows must be in time order: `fe_data_prep()` keeps each
  patient's rows in their input order, so sort by time within patient first.
- Direct standardization merges rows with equal linear predictors before
  summing, so it runs in O(providers x distinct linear predictors).
- `SR_output()` rates are on the 0-1 scale (pprof uses %), with `N_pop` as
  the denominator.
