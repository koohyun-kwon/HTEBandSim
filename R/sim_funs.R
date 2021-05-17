#' Noise generation
#'
#' @param n.sim number of simulation to use.
#' @inheritParams cvar.fun
#' @param dist distribution for noise function; either normal (\code{"norm"}) or log-normal distribution (\code{"lnorm"}).
#' @param scale scale parameter, a positive constant multiplied to the noise.
#'
#' @return matrix of noise terms with \code{nrow = length(x)} and \code{ncol = n.sim}.
#' @export
#'
#' @examples
#' x <- seq(-1, 1, length.out = 100)
#' eps_gen(150, "het.AK", "lnorm", x, 1/2)
eps_gen <- function(n.sim, cvar.spec, dist = c("norm", "lnorm"), x, scale){

  dist <- match.arg(dist)
  len.all <- n.sim * length(x)
  cvar.all <- rep(cvar.fun(cvar.spec, x), n.sim)
  rn <-
    if(dist == "norm"){
      stats::rnorm(len.all)
    }else if(dist == "lnorm"){
      (stats::rlnorm(len.all) - exp(1/2)) / sqrt(exp(2) - exp(1)) # location-scale normalization to yield 0 mean and unit variance
    }

  res <- matrix(sqrt(cvar.all) * rn * scale, length(x), n.sim)
  return(res)
}

#' True regression function value
#'
#' @param reg.name class of regression functions to be used; \code{"AK"} denotes regression functions
#' in Section 5 of Armstrong and Kolesár (2020); \code{"AK.bin"} denotes denotes the case where
#' \eqn{y = 1(f(x) + u > 0)}, using regression functions
#' in Section 5 of Armstrong and Kolesár (2020), where \eqn{u} is a normal random variable.
#' @param reg.spec a number denoting the type of regression function among the class of functions chosen
#' in \code{reg.name}.
#' @inheritParams f_AK
#' @inheritParams eps_gen
#'
#' @return a vector of true regression function values evaluated at \code{x}
#' @export
#'
#' @examples
#' x <- seq(-1, 1, length.out = 50)
#' true_f("AK", 1, 1, x)
#' true_f("AK.bin",1, 1, x, "het.AK", 1/2)
true_f <- function(reg.name = c("AK", "AK.bin"), reg.spec, M, x, cvar.spec = NULL, scale = NULL){

    f.x <-
      if(reg.name == "AK"){
        f_AK(reg.spec, x, M)
      }else if(reg.name == "AK.bin"){
        stats::pnorm(f_AK(reg.spec, x, M) / (sqrt(cvar.fun(cvar.spec, x)) * scale))  # Normal distribution is assumed
      }

    return(f.x)
}

#' Outcome generation
#'
#' @inheritParams eps_gen
#' @inheritParams true_f
#'
#' @return matrix of outcome variables with \code{nrow = length(x)} and \code{ncol = n.sim}.
#' @export
#'
#' @examples
#' x <- seq(-1, 1, length.out = 50)
#' y_gen(70, "AK", 3, 1, x, "het.AK", "lnorm", 1/2)
#' y_gen(70, "AK.bin", 2, 1, x, "het.AK", "norm", 1/2)
y_gen <- function(n.sim, reg.name = c("AK", "AK.bin"), reg.spec, M, x, cvar.spec, dist, scale){

  reg.name <- match.arg(reg.name)
  if(reg.name == "AK.bin") dist <- "norm"  # Only normal errors are supported for binary outcome case

  f.x <- true_f(reg.name, reg.spec, M, x, cvar.spec, scale)
  eps <- eps_gen(n.sim, cvar.spec, dist, x, scale)

  y.mat <- matrix(rep(f.x, n.sim), length(x), n.sim) + eps
  if(reg.name == "AK.bin"){
    y.mat <- (y.mat > 0)^2  # Transforming logicals into numeric values
  }

  return(y.mat)
}

#' Observation generation
#'
#' @param n number of observations
#' @inheritParams y_gen
#' @param x.spec distribution of the regressor, either uniform distribution (\code{"unif"}) or
#' normalized beta(2, 2) distribution (\code{"beta"}), both with the support of \eqn{[-1, 1]}.
#' @param x.l lower end of the support of the regressor.
#' @param x.u upper end of the support of the regressor.
#'
#' @return list of following components
#' \describe{
#' \item{x}{vector of regressors}
#' \item{y}{matrix of outcome variables with \code{nrow = n} and \code{ncol = n.sim}}
#' }
#' @export
#'
#' @examples obs_gen(100, 100, "AK.bin", 2, 1, "het.AK", "norm", 1/2, "beta")
obs_gen <- function(n, n.sim, reg.name = c("AK", "AK.bin"), reg.spec, M, cvar.spec, dist, scale,
                    x.spec = c("unif", "beta"), x.l = -1, x.u = 1){

  x.spec <- match.arg(x.spec)
  x <- if(x.spec == "unif"){
    stats::runif(n, min = x.l, max = x.u)
  }else if(x.spec == "beta"){
   (x.u - x.l) * stats::rbeta(n, 2, 2) + x.l
  }

  y.mat <- y_gen(n.sim, reg.name, reg.spec, M, x, cvar.spec, dist, scale)

  res <- list(x = x, y = y.mat)
  return(res)
}

#' Coverage indicator generation
#'
#' @param eval a vector of evaluation points in the support of the regressor
#' @param cb.l a vector of lower confidence band values
#' @param cb.u a vector of upper confidence band values
#' @inheritParams true_f
#'
#' @return a value of either 1 or 0, indicating coverage or non-coverage, respectively.
#' @export
covind_gen <- function(eval, cb.l, cb.u, reg.name = c("AK", "AK.bin"), reg.spec, M, cvar.spec = NULL, scale = NULL){

  f.x <- true_f(reg.name, reg.spec, M, eval, cvar.spec, scale)
  cov.ind.l <- cb.l <= f.x
  cov.ind.u <- cb.u >= f.x

  noncov.ind <- (FALSE %in% cov.ind.l) | (FALSE %in% cov.ind.u)
  cov.ind <- as.numeric(!noncov.ind)

  return(cov.ind)
}


#' Weighted length of a confidence band
#'
#' Computes weighted length of a confidence band, where the weights are given
#' by the design of the regressor.
#'
#' @inheritParams covind_gen
#' @param x.spec distribution of the regressor, which determines the weights. Supports
#' the uniform distribution (\code{"unif"}) and the normalized beta(2, 2) distribution (\code{"beta"}),
#' both with the support of \eqn{[-1, 1]}.
#'
#' @return a scalar value of weighted length.
#' @export
weighted_len <- function(eval, cb.l, cb.u, x.spec = c("unif", "beta")){

  x.spec <- match.arg(x.spec)

  cb.len <- cb.u - cb.l
  weights <-
    if(x.spec == "unif"){
      rep(1 / length(eval), length(eval))
    }else if(x.spec == "beta"){
      stats::dbeta(eval, 2, 2)
    }

  res <- sum(cb.len * weights) / sum(weights)
  return(res)
}

#' Coverage and length results generation
#'
#' Generates 1) coverage indicator, 2) weighted length, and 3) supremum length
#' given a confidence band.
#'
#' @inheritParams covind_gen
#' @inheritParams weighted_len
#'
#' @return a triplet of 1) coverage indicator, 2) weighted length, and 3) supremum length
#' @export
all_gen <- function(eval, cb.l, cb.u, reg.name = c("AK", "AK.bin"), reg.spec, M, cvar.spec = NULL, scale = NULL,
                    x.spec = c("unif", "beta")){

  cov.ind <- covind_gen(eval, cb.l, cb.u, reg.name, reg.spec, M, cvar.spec, scale)
  len.w <- weighted_len(eval, cb.l, cb.u, x.spec)
  len.sup <- max(cb.u - cb.l)

  return(c(cov.ind, len.w, len.sup))
}
