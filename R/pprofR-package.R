#' pprofR: shared modelling, testing and standardized measures for KECC quality measures
#'
#' One engine for the fixed provider effects logistic model used by the AOH,
#' PPPW, SFR and SRR measures. The workflow mirrors the measure code:
#'
#' 1. [fe_data_prep()] builds the `data_prep` object (binary or binomial
#'    outcome, optional patient identifier for clustered inference);
#' 2. [logis_fe()] fits the model (SerBIN or BAN);
#' 3. [test_fe()] tests provider effects, [summary_fe_covar()] tests covariate
#'    effects, [confint.logis_fe()] gives confidence intervals;
#' 4. [SR_output()] and [direct_exp()] compute indirect / direct standardized
#'    ratios and rates, and [test_SM()] tests a standardized measure against a
#'    reference value.
#'
#' Empirical-null calibration is deliberately not part of the package: the
#' measure code estimates the null location / scale and passes them to
#' [test_SM()].
#'
#' @keywords internal
#' @useDynLib pprofR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats coef confint vcov
"_PACKAGE"
