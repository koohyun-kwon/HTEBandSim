#' Debiased local polynomial estimation.
#'
#' A code provided by Gang Cheng at University of Washington (gangc at uw.edu)
#
#' @param X univariate X random variable.
#' @param Y univariate Y random variable.
#' @param h bandwidth for estimating the regression function.
#' @param tau ratio for bandwidths of estimating the derivative function over the bandwidth
#          for estimating the regression function.
#' @param x_seq x data point for evaluting the local polynomial estimator.
#'
#' @return a list with estimates with original local polynomial estimation (\code{lpr}),
#' and estimates with debiased local polynomial estimation (\code{lpr_de}).
#' @export
#'
#' @examples x <- stats::runif(500, min = -1, max = 1)
#' y <- x + rnorm(500, 0, 1/4)
#' eval <- seq(from = -0.9, to = 0.9, length.out = 20)
#' de_lpr(x, y, 0.1, x_seq = eval)
de_lpr <- function(X, Y, h, tau = 1, x_seq) {

  lpr = locpol::locPolSmootherC(X, Y, x_seq, h, deg  = 1,
                                kernel = locpol::gaussK)$beta0
  lpr_deriv = locpol::locPolSmootherC(X, Y, x_seq, h/tau, deg = 3,
                                      kernel = locpol::gaussK)$beta2
  lpr_de = lpr - 0.5 * h^2 * lpr_deriv
  return(list("lpr" = lpr, "lpr_de" = lpr_de))
}


#' Confidence Band Using Bootstrap
#'
#' Genereates a confidence band using the procedure described in
#' Cheng and Chen (2019).
#'
#' \code{B = 1000} is used in Cheng and Chen (2019).
#'
#' @param y dependent variable.
#' @param x independent variable.
#' @param eval vector of evaluation points.
#' @param B number of bootstrap simulations; default is \code{B = 1000}.
#' @param level confidence level used for confidence band; default is
#' \code{level = 0.95}.
#' @param fixed calculates the fixed length if \code{TRUE}; otherwise, CB is calculated
#' according to Remark 2 of Cheng and Chen (2019).
#'
#' @return a data frame containing index set and corresponding confidence band values
#' @export
#'
#' @references{
#'
#' \cite{Cheng,  G., and Y.-C. Chen. 2019.
#' "Nonparametric  inference  via  bootstrapping  the
#' debiased estimator."
#' Electronic Journal of Statistics, 13 (1): 2194–2256.}
#' }
#'
#' @examples x <- stats::runif(500, min = -1, max = 1)
#' y <- x + rnorm(500, 0, 1/4)
#' eval <- seq(from = -1, to = 1, length.out = 20)
#' CB_RBC(y, x, eval, 2)
#' CB_RBC(y, x, eval, 2, fixed = FALSE)
CB_RBC <- function(y, x, eval, B = 1000, level = 0.95, fixed = TRUE){

  est_res <- nprobust::lprobust(y, x, eval, bwselect = "mse-dpi")
  est_vec <- est_res$Estimate[, "tau.bc"]
  h.bw <- est_res$Estimate[, "h"]
  b.bw <- est_res$Estimate[, "b"]

  n <- length(y)

  max_d_vec <- numeric(B)

  if(fixed == TRUE){

    for(b in 1:B){

      b_ind <- sample(1:n, n, TRUE)
      y_b <- y[b_ind]
      x_b <- x[b_ind]

      est_res_b <- nprobust::lprobust(y_b, x_b, eval, h = h.bw, b = b.bw)
      est_vec_b <- est_res_b$Estimate[, "tau.bc"]

      max_d_vec[b] <- max(abs(est_vec_b - est_vec))
    }

    shat <- stats::quantile(max_d_vec, level)

    res_l <- est_vec - shat
    res_u <- est_vec + shat

  }else{

    for(b in 1:B){

      b_ind <- sample(1:n, n, TRUE)
      y_b <- y[b_ind]
      x_b <- x[b_ind]

      est_res_b <- nprobust::lprobust(y_b, x_b, eval, h = h.bw, b = b.bw)
      est_vec_b <- est_res_b$Estimate[, "tau.bc"]
      se_vec_b <- est_res_b$Estimate[, "se.rb"]

      max_d_vec[b] <- max(abs(est_vec_b - est_vec) / se_vec_b )
    }

    shat <- stats::quantile(max_d_vec, level)

    sd_vec <- est_res$Estimate[, "se.rb"]

    res_l <- est_vec - shat * sd_vec
    res_u <- est_vec + shat * sd_vec

  }

  res <- data.frame(eval = eval, cb.lower = res_l, cb.upper = res_u, h = h.bw)
  return(res)
}
