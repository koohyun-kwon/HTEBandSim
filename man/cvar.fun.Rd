% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DGP.R
\name{cvar.fun}
\alias{cvar.fun}
\title{Conditional variance function}
\usage{
cvar.fun(cvar.spec = c("homo", "het.AK"), x)
}
\arguments{
\item{cvar.spec}{specification for conditional variance function, one of \code{c("homo", "het.AK")}.}

\item{x}{a vector of regressors}
}
\value{
a vector of conditional variance values
}
\description{
Conditional variance function
}
\examples{
x <- seq(-1, 1, length.out = 100)
cvar.fun("homo", x)
cvar.fun("het.AK", x)
}
