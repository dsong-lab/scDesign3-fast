test_that("extract_para reuses mgcv prediction design for shared smooths", {
  skip_if_not_installed("mgcv")
  skip_if_not_installed("SingleCellExperiment")

  set.seed(301)
  n <- 70L
  g <- 3L
  dat <- data.frame(pseudotime = sort(stats::runif(n)))
  rownames(dat) <- paste0("cell", seq_len(n))
  eta <- 0.5 + cos(2 * pi * dat$pseudotime)
  count_mat <- matrix(
    stats::rnbinom(n * g, size = 4, mu = rep(exp(eta), g)),
    nrow = n,
    ncol = g
  )
  rownames(count_mat) <- rownames(dat)
  colnames(count_mat) <- paste0("gene", seq_len(g))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = t(count_mat))
  )

  fit <- fit_marginal(
    data = list(count_mat = count_mat, dat = dat, filtered_gene = NULL),
    mu_formula = "s(pseudotime, bs = 'cr', k = 8)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    use_scglm = "never"
  )

  new_covariate <- data.frame(pseudotime = seq(0.02, 0.98, length.out = 50))
  rownames(new_covariate) <- paste0("newcell", seq_len(nrow(new_covariate)))
  fit_list <- lapply(fit, function(x) x$fit)
  removed_cell_list <- lapply(fit, function(x) x$removed_cell)
  cache <- .scdesign3_prepare_predict_cache(
    marginal_list = fit_list,
    qc_gene_idx = seq_len(g),
    new_covariate = new_covariate,
    removed_cell_list = removed_cell_list
  )

  para <- extract_para(
    sce = sce,
    marginal_list = fit,
    n_cores = 1,
    family_use = "nb",
    new_covariate = new_covariate,
    data = dat
  )

  direct_mean <- sapply(fit, function(x) {
    stats::predict(x$fit, type = "response", newdata = new_covariate)
  })
  direct_theta <- sapply(fit, function(x) {
    rep(1 / x$fit$family$getTheta(TRUE), nrow(new_covariate))
  })

  expect_false(is.null(cache))
  expect_equal(unname(para$mean_mat), unname(direct_mean), tolerance = 1e-8)
  expect_equal(unname(para$sigma_mat), unname(direct_theta), tolerance = 1e-8)
})
