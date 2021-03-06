% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival-aareg-tidiers.R
\name{tidy.aareg}
\alias{tidy.aareg}
\alias{aareg_tidiers}
\title{Tidy a(n) aareg object}
\usage{
\method{tidy}{aareg}(x, ...)
}
\arguments{
\item{x}{An \code{aareg} object returned from \code{\link[survival:aareg]{survival::aareg()}}.}

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
A \link[tibble:tibble]{tibble::tibble} with one row for each coefficient and columns:

\item{term}{name of coefficient}
\item{estimate}{estimate of the slope}
\item{statistic}{test statistic for coefficient}
\item{std.error}{standard error of statistic}
\item{robust.se}{robust version of standard error estimate (only when
\code{x} was called with \code{dfbeta = TRUE})}
\item{z}{z score}
\item{p.value}{p-value}
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

library(survival)

afit <- aareg(
  Surv(time, status) ~ age + sex + ph.ecog,
  data = lung,
  dfbeta = TRUE
)

tidy(afit) 

}
\seealso{
\code{\link[=tidy]{tidy()}}, \code{\link[survival:aareg]{survival::aareg()}}

Other aareg tidiers: \code{\link{glance.aareg}}

Other survival tidiers: \code{\link{augment.coxph}},
  \code{\link{augment.survreg}},
  \code{\link{glance.aareg}}, \code{\link{glance.cch}},
  \code{\link{glance.coxph}}, \code{\link{glance.pyears}},
  \code{\link{glance.survdiff}},
  \code{\link{glance.survexp}},
  \code{\link{glance.survfit}},
  \code{\link{glance.survreg}}, \code{\link{tidy.cch}},
  \code{\link{tidy.coxph}}, \code{\link{tidy.pyears}},
  \code{\link{tidy.survdiff}}, \code{\link{tidy.survexp}},
  \code{\link{tidy.survfit}}, \code{\link{tidy.survreg}}
}
\concept{aareg tidiers}
\concept{survival tidiers}
