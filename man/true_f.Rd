% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_funs.R
\name{true_f}
\alias{true_f}
\title{True regression function value}
\usage{
true_f(
  reg.name = c("AK", "AK.bin"),
  reg.spec,
  M,
  x,
  cvar.spec = NULL,
  scale = NULL
)
}
\arguments{
\item{reg.name}{class of regression functions to be used; \code{"AK"} denotes regression functions
in Section 5 of Armstrong and Kolesár (2020); \code{"AK.bin"} denotes denotes the case where
\eqn{y = 1(f(x) + u > 0)}, using regression functions
in Section 5 of Armstrong and Kolesár (2020), where \eqn{u} is a normal random variable.}

\item{reg.spec}{a number denoting the type of regression function among the class of functions chosen
in \code{reg.name}.}

\item{M}{bound on the second derivative}

\item{x}{a vector of x values.}

\item{cvar.spec}{specification for conditional variance function, one of \code{c("homo", "het.AK")}.}

\item{scale}{scale parameter, a positive constant multiplied to the noise.}
}
\value{
a vector of true regression function values evaluated at \code{x}
}
\description{
True regression function value
}
\examples{
x <- seq(-1, 1, length.out = 50)
true_f("AK", 1, 1, x)
true_f("AK.bin",1, 1, x, "het.AK", 1/2)
}
