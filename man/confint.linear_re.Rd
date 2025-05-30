% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confint.linear_re.R
\name{confint.linear_re}
\alias{confint.linear_re}
\title{Get confidence intervals for provider effects or standardized measures from a fitted \code{linear_re} object}
\usage{
\method{confint}{linear_re}(
  object,
  parm,
  level = 0.95,
  option = "SM",
  stdz = "indirect",
  alternative = "two.sided",
  ...
)
}
\arguments{
\item{object}{a model fitted from \code{linear_re}.}

\item{parm}{specify a subset of providers for which confidence intervals are given.
By default, all providers are included. The class of \code{parm} should match the class of the provider IDs.}

\item{level}{the confidence level. The default value is 0.95.}

\item{option}{a character string specifying whether the confidence intervals
should be provided for provider effects or standardized measures:
\itemize{
\item {\code{"alpha"}} provider effect.
\item {\code{"SM"}} standardized measures.
}}

\item{stdz}{a character string or a vector specifying the standardization method
if \code{option} includes \code{"SM"}. See \code{stdz} argument in \code{\link{SM_output.linear_re}}.}

\item{alternative}{a character string specifying the alternative hypothesis, must be one of
\code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.
Note that \code{"alpha"} for argument \code{option} only supports \code{"two.sided"}.}

\item{\dots}{additional arguments that can be passed to the function.}
}
\value{
A list of data frames containing the confidence intervals based on the values of \code{option} and \code{stdz}.
\item{CI.alpha}{Confidence intervals for provider effects if \code{option} includes \code{"alpha"}.}
\item{CI.indirect}{Confidence intervals for indirect standardized differences if \code{option} includes \code{"SM"} and \code{stdz} includes \code{"indirect"}.}
\item{CI.direct}{Confidence intervals for direct standardized differences if \code{option} includes \code{"SM"} and \code{stdz} includes \code{"direct"}.}
}
\description{
Provide confidence intervals for provider effects or standardized measures from a random effect linear model.
}
\examples{
data(ExampleDataLinear)
outcome <- ExampleDataLinear$Y
ProvID <- ExampleDataLinear$ProvID
covar <- ExampleDataLinear$Z
fit_re <- linear_re(Y = outcome, Z = covar, ProvID = ProvID)
confint(fit_re)

}
