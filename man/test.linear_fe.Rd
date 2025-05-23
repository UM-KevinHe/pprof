% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.linear_fe.R
\name{test.linear_fe}
\alias{test.linear_fe}
\title{Conduct hypothesis testing for provider effects from a fitted \code{linear_fe} object}
\usage{
\method{test}{linear_fe}(fit, parm, level = 0.95, null = "median", alternative = "two.sided", ...)
}
\arguments{
\item{fit}{a model fitted from \code{linear_fe}.}

\item{parm}{specifies a subset of providers for which confidence intervals are to be given.
By default, all providers are included. The class of \code{parm} should match the class of the provider IDs.}

\item{level}{the confidence level during the hypothesis test, meaning a significance level of \eqn{1 - \text{level}}.
The default value is 0.95.}

\item{null}{a character string or a number defining the null hypothesis for the provider effects.
The default value is \code{"median"}. The possible values are:
\itemize{
\item{\code{"median"}}: The median of the provider effect estimates (\eqn{\hat{\gamma}_i}).
\item{\code{"mean"}}: The weighted average of the provider effect estimates (\eqn{\hat{\gamma}_i}), where the weights correspond to the sample size of each provider.
\item{numeric}: A user-defined numeric value representing the null hypothesis.
}}

\item{alternative}{a character string specifying the alternative hypothesis, must be one of
\code{"two.sided"} (default), \code{"greater"}, or \code{"less"}.}

\item{\dots}{additional arguments that can be passed to the function.}
}
\value{
A data frame containing the results of the hypothesis test, with the following columns:
\item{flag}{a flagging indicator where \code{1} means statistically higher than expected
and \code{-1} means statistically lower than expected.}
\item{p-value}{the p-value of the hypothesis test.}
\item{stat}{the test statistic.}
\item{Std.Error}{the standard error of the provider effect estimate.}
}
\description{
Conduct hypothesis tests on provider effects and identify outlying providers for a fixed effect linear model.
}
\details{
The function identifies outlying providers based on hypothesis test results.
For two-sided tests, \code{1} indicates performance significantly higher than expected, \code{-1} indicates lower,
For one-sided tests, \code{1} (right-tailed) or \code{-1} (left-tailed) flags are used.
Providers whose performance falls within the central range are flagged as \code{0}.
Outlying providers are determined by the test statistic falling beyond the threshold based on the significance level \eqn{1 - \text{level}}.
}
\examples{
data(ExampleDataLinear)
outcome <- ExampleDataLinear$Y
covar <- ExampleDataLinear$Z
ProvID <- ExampleDataLinear$ProvID
fit_linear <- linear_fe(Y = outcome, Z = covar, ProvID = ProvID)
test(fit_linear)

}
