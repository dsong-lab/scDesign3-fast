test_that("cached Gaussian copula scores match dense mvtnorm likelihood", {
  skip_if_not_installed("mvtnorm")

  set.seed(201)
  n <- 120L
  p <- 12L
  z <- matrix(stats::rnorm(n * p), nrow = n)
  colnames(z) <- paste0("gene", seq_len(p))

  cor_mat <- cal_cor(
    norm.mat = z,
    important_feature = rep(TRUE, p),
    score_mat = z
  )
  cache <- attr(cor_mat, "scdesign3_gaussian_factor", exact = TRUE)
  expect_true(!is.null(cache))

  old_loglik <- sum(mvtnorm::dmvnorm(
    x = z,
    mean = rep(0, p),
    sigma = as.matrix(cor_mat),
    log = TRUE
  )) - sum(rowSums(stats::dnorm(z, log = TRUE)))
  npar <- (as.integer(sum(as.matrix(cor_mat) != 0)) - p) / 2

  expect_equal(attr(cor_mat, "scdesign3_gaussian_loglik", exact = TRUE), old_loglik, tolerance = 1e-8)
  expect_equal(cal_aic(z, cor_mat, ind = FALSE), -2 * old_loglik + 2 * npar, tolerance = 1e-8)
  expect_equal(cal_bic(z, cor_mat, ind = FALSE), -2 * old_loglik + log(n) * npar, tolerance = 1e-8)
})

test_that("cached Gaussian copula sampling uses an exact reusable factor", {
  set.seed(202)
  n <- 20L
  p <- 8L
  z <- matrix(stats::rnorm(100L * p), nrow = 100L)
  colnames(z) <- paste0("gene", seq_len(p))
  cor_cached <- cal_cor(z, rep(TRUE, p), score_mat = z)
  cache <- attr(cor_cached, "scdesign3_gaussian_factor", exact = TRUE)

  expect_true(!is.null(cache))
  expect_equal(crossprod(cache$root), as.matrix(cor_cached), tolerance = 1e-10, ignore_attr = TRUE)
  set.seed(203)
  new <- sampleMVN(n = n, Sigma = cor_cached, n_cores = 1L, fastmvn = FALSE)

  expect_equal(dim(new), c(n, p))
  expect_true(all(new > 0 & new < 1))
})

test_that("Gaussian copula detects negatively correlated genes as correlated", {
  sigma <- matrix(c(1, -0.4, 0, -0.4, 1, 0, 0, 0, 1), nrow = 3L)
  colnames(sigma) <- rownames(sigma) <- c("a", "b", "c")

  expect_equal(.scdesign3_correlated_gene_names(sigma), c("a", "b"))
})

test_that("block-factor Gaussian copula scores and samples", {
  set.seed(206)
  n <- 90L
  p <- 12L
  f <- matrix(stats::rnorm(n * 2L), nrow = n)
  load <- matrix(stats::rnorm(p * 2L, sd = 0.4), nrow = p)
  z <- f %*% t(load) + matrix(stats::rnorm(n * p, sd = 0.7), nrow = n)
  colnames(z) <- paste0("gene", seq_len(p))

  fit <- .scdesign3_fit_gaussian_copula(
    norm.mat = z,
    important_feature = rep(TRUE, p),
    gaussian_copula = "block_factor",
    gaussian_copula_rank = 2L,
    gaussian_copula_block_size = 4L,
    gaussian_copula_block_shrinkage = 0.05,
    gaussian_copula_sketch_size = 40L,
    score_mat = z
  )

  expect_s3_class(fit, "scdesign3_gaussian_fast_copula")
  expect_equal(fit$method, "block_factor")
  expect_true(length(fit$blocks) >= 3L)
  expect_true(max(vapply(fit$blocks, length, integer(1))) <= 4L)
  expect_true(is.finite(cal_aic(z, fit, ind = FALSE)))
  expect_true(is.finite(cal_bic(z, fit, ind = FALSE)))

  set.seed(207)
  q <- sampleMVN(n = 20L, Sigma = fit, n_cores = 1L, fastmvn = FALSE)
  expect_equal(dim(q), c(20L, p))
  expect_equal(colnames(q), colnames(z))
  expect_true(all(q > 0 & q < 1))
})
