test_that("cached mgcv spline setup matches direct per-gene mgcv fits", {
  skip_if_not_installed("mgcv")

  set.seed(201)
  n <- 80L
  g <- 4L
  dat <- data.frame(pseudotime = stats::runif(n))
  rownames(dat) <- paste0("cell", seq_len(n))
  eta <- 1 + sin(2 * pi * dat$pseudotime)
  count_mat <- matrix(
    stats::rnbinom(n * g, size = 3, mu = rep(exp(eta), g)),
    nrow = n,
    ncol = g
  )
  rownames(count_mat) <- rownames(dat)
  colnames(count_mat) <- paste0("gene", seq_len(g))

  cached <- fit_marginal(
    data = list(count_mat = count_mat, dat = dat, filtered_gene = NULL),
    mu_formula = "s(pseudotime, bs = 'cr', k = 8)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    use_scglm = "never"
  )

  direct <- lapply(seq_len(g), function(j) {
    dat_j <- dat
    dat_j$gene <- count_mat[, j]
    mgcv::gam(gene ~ s(pseudotime, bs = "cr", k = 8), data = dat_j, family = "nb", method = "REML")
  })

  cached_pred <- unlist(lapply(cached, function(z) as.numeric(stats::predict(z$fit, type = "response"))))
  direct_pred <- unlist(lapply(direct, function(z) as.numeric(stats::predict(z, type = "response"))))
  cached_theta <- vapply(cached, function(z) z$fit$family$getTheta(TRUE), numeric(1))
  direct_theta <- vapply(direct, function(z) z$family$getTheta(TRUE), numeric(1))

  expect_equal(unname(cached_pred), unname(direct_pred), tolerance = 1e-7)
  expect_equal(unname(cached_theta), unname(direct_theta), tolerance = 1e-7)
})

test_that("gamlss ga prefit cache reuses shared spline setup", {
  skip_if_not_installed("mgcv")

  set.seed(202)
  n <- 60L
  dat <- data.frame(
    pseudotime = stats::runif(n),
    Y.var = 0
  )
  form <- stats::as.formula("Y.var ~ s(pseudotime, bs = 'cr', k = 6)")
  control <- ga.control(method = "REML")
  cache_env <- new.env(parent = emptyenv())
  old <- options(.scdesign3_gamlss_cache_env = cache_env)
  on.exit(options(old), add = TRUE)

  G1 <- .scdesign3_gamlss_ga_prefit(form, dat, control)
  G2 <- .scdesign3_gamlss_ga_prefit(form, dat, control)

  expect_length(ls(cache_env), 1L)
  expect_identical(G1, G2)
})
