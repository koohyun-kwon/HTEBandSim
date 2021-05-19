test_that("coverage indicator and length functions work", {
  # obs <- obs_gen(n = 500, n.sim = 100, "AK", reg.spec = 2, M = 1, "het.AK", "norm", scale = 1/2, "beta")
  # x <- obs$x
  # y <- obs$y
  # y <- y[, 1]
  #
  # cb.res <- HTEBand::NpregBand(y, x, 1, 0.95, "H", n.eval = 50, q.int = 0.05)
  # covind_gen(cb.res$eval, cb.res$cb.lower, cb.res$cb.upper, "AK", reg.spec = 2, M = 1)
  # weighted_len(cb.res$eval, cb.res$cb.lower, cb.res$cb.upper, x.spec = "beta")
  # all_gen(cb.res$eval, cb.res$cb.lower, cb.res$cb.upper, cb.res$h.t, cb.res$h.c, "AK", reg.spec = 2, M = 1, x.spec = "beta")
})


test_that("coverage indicator and length functions work:RBC", {
  # obs <- obs_gen(n = 500, n.sim = 100, "AK", reg.spec = 2, M = 1, "het.AK", "norm", scale = 1/2, "beta")
  # x <- obs$x
  # y <- obs$y
  # y <- y[, 1]
  # eval <- seq(from = -1, to = 1, length.out = 20)
  #
  # cb.res <- CB_RBC(y, x, eval, B = 2, fixed = FALSE, bwselect = "ce-dpi")
  # cb.res2 <- CB_RBC(y, x, eval, B = 2, fixed = FALSE, bwselect = "mse-dpi")
  # all_gen(cb.res$eval, cb.res$cb.lower, cb.res$cb.upper, cb.res$h, cb.res$h, "AK", reg.spec = 2, M = 1, x.spec = "beta")
})
