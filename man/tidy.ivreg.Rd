% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aer-tidiers.R
\name{tidy.ivreg}
\alias{tidy.ivreg}
\alias{ivreg_tidiers}
\alias{aer_tidiers}
\title{Tidy a(n) ivreg object}
\usage{
\method{tidy}{ivreg}(x, conf.int = FALSE, conf.level = 0.95,
  exponentiate = FALSE, ...)
}
\arguments{
\item{x}{An \code{ivreg} object created by a call to \code{\link[AER:ivreg]{AER::ivreg()}}.}

\item{conf.int}{Logical indicating whether or not to include a confidence
interval in the tidied output. Defaults to \code{FALSE}.}

\item{conf.level}{The confidence level to use for the confidence interval
if \code{conf.int = TRUE}. Must be strictly greater than 0 and less than 1.
Defaults to 0.95, which corresponds to a 95 percent confidence interval.}

\item{exponentiate}{Logical indicating whether or not to exponentiate the
the coefficient estimates. This is typical for logistic and multinomial
regressions, but a bad idea if there is no log or logit link. Defaults
to \code{FALSE}.}

\item{...}{Additional arguments. Not used. Needed to match generic
signature only. \strong{Cautionary note:} Misspelled arguments will be
absorbed in \code{...}, where they will be ignored. If the misspelled
argument has a default value, the default value will be used.
For example, if you pass \code{conf.lvel = 0.9}, all computation will
proceed using \code{conf.level = 0.95}. Additionally, if you pass
\code{newdata = my_tibble} to an \code{\link[=augment]{augment()}} method that does not
accept a \code{newdata} argument, it will use the default value for
the \code{data} argument.}
}
\value{
A \code{\link[tibble:tibble]{tibble::tibble()}} with one row for each term in the
regression. The tibble has columns:

\item{term}{The name of the regression term.}
\item{estimate}{The estimated value of the regression term.}
\item{std.error}{The standard error of the regression term.}
\item{statistic}{The value of a statistic, almost always a T-statistic,
to use in a hypothesis that the regression term is non-zero.}
\item{p.value}{The two-sided p-value associated with the observed
statistic.}
\item{conf.low}{The low end of a confidence interval for the regression
term. Included only if \code{conf.int = TRUE}.}
\item{conf.high}{The high end of a confidence interval for the regression
term. Included only if \code{conf.int = TRUE}.}
}
\description{
Tidy summarizes information about the components of a model.
A model component might be a single term in a regression, a single
hypothesis, a cluster, or a class. Exactly what tidy considers to be a
model component varies cross models but is usually self-evident.
If a model has several distinct types of components, you will need to
specify which components to return.
}
\examples{

library(AER)

data("CigarettesSW", package = "AER")
ivr <- ivreg(
  log(packs) ~ income | population,
  data = CigarettesSW,
  subset = year == "1995"
)

summary(ivr)

tidy(ivr)
tidy(ivr, conf.int = TRUE)
tidy(ivr, conf.int = TRUE, exponentiate = TRUE)

augment(ivr)

glance(ivr)

}
\seealso{
\code{\link[=tidy]{tidy()}}, \code{\link[AER:ivreg]{AER::ivreg()}}

Other ivreg tidiers: \code{\link{augment.ivreg}},
  \code{\link{glance.ivreg}}
}
\concept{ivreg tidiers}
