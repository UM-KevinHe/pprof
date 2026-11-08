#' Standardized measures of a fitted provider model
#'
#' Direct and indirect standardized measures (ratios, rates or differences,
#' depending on the model) for every provider of a fitted risk-adjustment
#' model, with their standard errors where the model provides them.
#'
#' @param fit a fitted provider model.
#' @param ... arguments passed to the method.
#' @return A data frame with one row per provider; see the methods.
#' @export
standardize <- function(fit, ...) {
  UseMethod("standardize")
}

#' Tests of the provider effects
#'
#' Tests the effect of every provider against a reference value and flags
#' the providers whose effect differs from it.
#'
#' @param fit a fitted provider model.
#' @param ... arguments passed to the method.
#' @return A data frame with one row per provider and a fixed set of columns;
#'   see the methods.
#' @export
test <- function(fit, ...) {
  UseMethod("test")
}
