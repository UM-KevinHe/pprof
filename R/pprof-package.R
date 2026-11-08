#' pprof: provider profiling by components
#'
#' Components for profiling health care providers, organised by analysis
#' step: data preparation ([pp_data()], [pp_check()]), risk-adjustment
#' models, standardized measures ([standardize()]), tests of provider effects
#' ([test_providers()]) and covariates ([test_covariates()]), reference
#' distributions, reliability, variable selection and plots. Every quality
#' measure is assembled from these components; nothing in the package is
#' specific to one measure.
#'
#' @section Conventions:
#' * Every function takes a data frame and column names (`y`, `x`,
#'   `provider`, `patient`, `n`; `time`, `start`, `event` for survival
#'   outcomes).
#' * Provider identifiers are handled as character strings and rows are
#'   sorted by provider, then patient, in the C locale; a patient's rows keep
#'   their input order, which the stabilized tests rely on.
#' * Rates are proportions on the 0-1 scale.
#' * Flags are `1` (above the reference), `-1` (below), `0` (not different)
#'   or `NA` (not tested); which tail is "worse" is a property of the measure
#'   and is left to the caller.
#'
#' @keywords internal
#' @useDynLib pprof, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef predict vcov
"_PACKAGE"
