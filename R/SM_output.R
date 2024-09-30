#' Generic Function for Calculating Standardized Measures
#'
#' `SM_output` is an S3 generic function designed to calculate standardized
#' measures. It dispatches to the appropriate method based on the class of the
#' input (`fit`), ensuring the correct method is applied for different types of models.
#'
#' @param fit the input object, typically a fitted model, for which standardized measures
#' are to be calculated. The method applied depends on the class of this object.
#' @param ... Additional arguments that can be passed to specific methods.
#'
#' @return the return varies depending on the method implemented for the
#' class of the input object.

SM_output <- function(fit, ...) {
  UseMethod("SM_output")
}
