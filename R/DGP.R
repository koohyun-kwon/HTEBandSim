#' Function \eqn{s(x)} in Armstrong and Kolesar (2020)
#'
#' Function \eqn{s(x)} as defined in Section 5 of Armstrong and Kolesár (2020).
#'
#' @param x a vector of x values
#'
#' @export
#'
#' @references{
#'
#' \cite{Armstrong, Timothy B., and Michal Kolesár. 2020.
#' "Simple and Honest Confidence Intervals in Nonparametric Regression."
#' Quantitative Economics 11 (1): 1–39.}
#'
#' }
#'
#' @examples
#' x <- runif(100, -1, 1)
#' s_fun(x)
s_fun <- function(x){

  res <- (pmax(x, 0))^2

  return(res)
}

#' Regression functions in Armstrong and Kolesar (2020)
#'
#' Functions \eqn{f_1(x)}, \eqn{f_2(x)}, and \eqn{f_3(x)} as defined in Section 5 of Armstrong and Kolesár (2020).
#'
#' @param spec which specification to use; choose among {1,2,3}.
#' @param x a vector of x values.
#' @param M bound on the second derivative
#'
#' @export
#'
#' @references{
#'
#' \cite{Armstrong, Timothy B., and Michal Kolesár. 2020.
#' "Simple and Honest Confidence Intervals in Nonparametric Regression."
#' Quantitative Economics 11 (1): 1–39.}
#'
#' }
#'
#' @examples x <- runif(100, -1, 1)
#' M <- 2
#' f_AK(1, x, M)
#' f_AK(2, x, M)
#' f_AK(3, x, M)
f_AK <- function(spec, x, M){

  res <-
    if(spec == 1){
      (M / 2) * (x^2 - 2 * s_fun(abs(x) - 0.25))
    }else if(spec == 2){
      (M / 2) * (x^2 - 2 * s_fun(abs(x) - 0.2) + 2 * s_fun(abs(x) - 0.5) - 2 * s_fun(abs(x) - 0.65))
    }else if(spec == 3){
      (M / 2) * ((x + 1)^2 - 2 * s_fun(x + 0.2) + 2 * s_fun(x - 0.2) - 2 * s_fun(x - 0.4)
                 + 2 * s_fun(x - 0.7) - 0.92)
    }

  return(res)
}



#' Berry, Caroll, and Ruppert (2002) Regression Function
#'
#' @param x a vector of x values.
#' @param M an additional parameter.
#'
#' @export
#'
#' @references{
#'
#' \cite{S. M. Berry, R. J. Carroll, and D. Ruppert. 2002.
#' "Bayesian smoothing and regression splines for measurement error problems."
#' Journal of the American Statistical Association 97 (457): 160–169.}
#'
#' }
#'
#' @examples x <- runif(100, -1, 1)
#' f.BCR(x)
f.BCR <- function(x, M = 3 * pi / 2){

  res <- sin(M * x) / (1 + 18 * x^2 * (sign(x) + 1))

  return(res)
}


#' Conditional variance function
#'
#' @param cvar.spec specification for conditional variance function, one of \code{c("homo", "het.AK")}.
#' @param x a vector of regressors
#'
#' @return a vector of conditional variance values
#' @export
#'
#' @examples
#' x <- seq(-1, 1, length.out = 100)
#' cvar.fun("homo", x)
#' cvar.fun("het.AK", x)
cvar.fun <- function(cvar.spec = c("homo", "het.AK"), x){

  cvar.spec <- match.arg(cvar.spec)

  n <- length(x)

  res <- if(cvar.spec == "homo"){

      rep(1, n)
  }else if(cvar.spec == "het.AK"){

      (1 + sqrt(abs(x)))^2  # as in Appendix F of Armstrong and Kolesár (2020, QE)
  }

  return(res)
}
