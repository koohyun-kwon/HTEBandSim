#' Noise generation
#'
#' @param n.sim number of simulation to use.
#' @inheritParams cvar.fun
#' @param dist distribution for noise function; choose among normal (\code{"norm"}),
#' log-normal distribution (\code{"lnorm"}) or logistic (\code{"logis"}).
#' @param scale scale parameter, a positive constant multiplied to the noise.
#'
#' @return matrix of noise terms with \code{nrow = length(x)} and \code{ncol = n.sim}.
#' @export
#'
#' @examples
#' x <- seq(-1, 1, length.out = 100)
#' eps_gen(150, "het.AK", "lnorm", x, 1/2)
eps_gen <- function(n.sim, cvar.spec, dist = c("norm", "lnorm", "logis"),
                    x, scale){

  dist <- match.arg(dist)
  len.all <- n.sim * length(x)
  cvar.all <- rep(cvar.fun(cvar.spec, x), n.sim)
  rn <-
    if(dist == "norm"){
      stats::rnorm(len.all)
    }else if(dist == "lnorm"){
      # My CB doesn't work well for this specification
      # Kurtosis = 76 >> 3 for normal => dist'n is too fat-tailed
      # In comparison, logistic noise has kurtosis = 4, moderately larger than 3

      # location-scale normalization to yield 0 mean and unit variance
      (stats::rlnorm(len.all) - exp(1/2)) / sqrt(exp(2) - exp(1))
    }else if(dist == "logis"){
      stats::rlogis(len.all, 0, sqrt(3) / pi)
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
#' @param custom.fx a vector of custom f(x) values
#'
#' @return a value of either 1 or 0, indicating coverage or non-coverage, respectively.
#' @export
covind_gen <- function(eval, cb.l, cb.u, reg.name = c("AK", "AK.bin", "custom"),
                       reg.spec, M, cvar.spec = NULL, scale = NULL,
                       custom.fx = NULL){

  if(reg.name != "custom"){
    f.x <- true_f(reg.name, reg.spec, M, eval, cvar.spec, scale)
  }else{
    f.x <- custom.fx
  }
  cov.ind.l <- cb.l <= f.x
  cov.ind.u <- cb.u >= f.x

  noncov.ind <- (FALSE %in% cov.ind.l) | (FALSE %in% cov.ind.u)
  cov.ind <- as.numeric(!noncov.ind)

  return(cov.ind)
}

#' Non-coverage "length" calculation
#'
#' @param eval a vector of evaluation points in the support of the regressor
#' @param cb.l a vector of lower confidence band values
#' @param cb.u a vector of upper confidence band values
#' @inheritParams true_f
#' @param custom.fx a vector of custom f(x) values
#'
#' @return a vector of length \code{2 * len(eval)}, with a vector of lower non-coverage length
#' and a vector of upper non-coverage length
#' @export
noncov_len <- function(eval, cb.l, cb.u, reg.name = c("AK", "AK.bin", "custom"),
                       reg.spec, M, cvar.spec = NULL, scale = NULL,
                       custom.fx = NULL){

  if(reg.name != "custom"){
    f.x <- true_f(reg.name, reg.spec, M, eval, cvar.spec, scale)
  }else{
    f.x <- custom.fx
  }

  cov.ind.l <- as.numeric(cb.l <= f.x)
  cov.ind.u <- as.numeric(cb.u >= f.x)

  noncov.len.l <- (cb.l - f.x) * (1 - cov.ind.l)
  noncov.len.u <- (f.x - cb.u) * (1 - cov.ind.u)

  return(c(noncov.len.l, noncov.len.u))
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
#' Generates 1) coverage indicator, 2) confidence band lengths at each evaluation points and 3) bandwidths used
#' at each evaluation points given a confidence band.
#'
#' @inheritParams covind_gen
#' @inheritParams weighted_len
#' @param h.t bandwidth used at each evaluation points for treated observations or for all observations
#' if there is no treated/control distinction.
#' @param h.c bandwidth used at each evaluation points for control observations;
#' it can be left unspecified if there is no treated/control distinction.
#'
#' @return a vector of 1) coverage indicator, 2) weighted length, 3) supremum length,
#' 4) a vector of confidence band lengths at each evaluation points
#' 5) a vector of bandwidths used at each evaluation points given a confidence band for treated and control observations,
#' 6) a vector of evaluation points
#' 7) a vector of "lower non-coverage length"
#' 8) a vector of "upper non-coverage length"
#' @export
all_gen <- function(eval, cb.l, cb.u, h.t, h.c = h.t,
                    reg.name = c("AK", "AK.bin", "custom"), reg.spec, M,
                    cvar.spec = NULL, scale = NULL, x.spec,
                    custom.fx = NULL){

  if(reg.name == "custom"){
    res.spec <- NULL
    M <- NULL
  }

  cov.ind <- covind_gen(eval, cb.l, cb.u, reg.name, reg.spec, M, cvar.spec, scale,
                        custom.fx)
  len.all <- cb.u - cb.l
  len.w <- weighted_len(eval, cb.l, cb.u, x.spec)
  len.sup <- max(len.all)
  len.noncov <- noncov_len(eval, cb.l, cb.u, reg.name, reg.spec, M, cvar.spec, scale,
                           custom.fx)

  res <- c(cov.ind, len.w, len.sup, len.all, h.t, h.c, eval, len.noncov)
  neval <- length(eval)
  res.name <- c("cov", "len.w", "len.sup", paste("len", 1:neval, sep = "."), paste("ht", 1:neval, sep = "."),
                paste("hc", 1:neval, sep = "."), paste("eval", 1:neval, sep = "."),
                paste("ncl.l", 1:neval, sep = "."), paste("ncl.u", 1:neval, sep = "."))
  names(res) <- res.name

  return(res)
}
