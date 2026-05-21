test_that("fit_marginal uses scGLM for shared categorical poisson models", {
  skip_if_not_installed("scGLM")
  skip_if_not_installed("SingleCellExperiment")

  set.seed(101)
  n <- 80L
  g <- 6L
  dat <- data.frame(group = factor(sample(c("a", "b"), n, TRUE)))
  rownames(dat) <- paste0("cell", seq_len(n))
  x <- stats::model.matrix(~ group, dat)
  beta <- matrix(stats::rnorm(ncol(x) * g, sd = 0.1), nrow = ncol(x))
  beta[1, ] <- log(3)
  count_mat <- matrix(stats::rpois(n * g, lambda = as.vector(exp(x %*% beta))), nrow = n)
  rownames(count_mat) <- rownames(dat)
  colnames(count_mat) <- paste0("gene", seq_len(g))

  fit <- fit_marginal(
    data = list(count_mat = count_mat, dat = dat, filtered_gene = NULL),
    mu_formula = "group",
    sigma_formula = "1",
    family_use = "poisson",
    n_cores = 2,
    use_scglm = "auto",
    scglm_batch_size = 2L
  )

  expect_s3_class(fit[[1]]$fit, "scdesign3_scglm")
  expect_identical(fit[[1]]$fit$backend, "scGLM::categorical_closed_form")
  expect_length(stats::predict(fit[[1]]$fit, type = "response"), n)
  expect_true(is.finite(stats::AIC(fit[[1]]$fit)))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = t(count_mat))
  )
  para <- extract_para(
    sce = sce,
    assay_use = "counts",
    marginal_list = fit,
    n_cores = 1,
    family_use = "poisson",
    new_covariate = dat,
    data = dat
  )
  expect_equal(unique(as.numeric(para$sigma_mat)), 1)
})

test_that("scglm_fit controls exact versus approximate automatic dispatch", {
  dat <- data.frame(
    group = factor(rep(c("a", "b"), each = 5L)),
    x = seq(0, 1, length.out = 10L)
  )

  expect_identical(
    .scdesign3_scglm_method("auto", "poisson", ~ group, dat, scglm_fit = "approximate"),
    "categorical_closed_form"
  )
  expect_identical(
    .scdesign3_scglm_method("auto", "poisson", ~ group, dat, scglm_fit = "exact"),
    "categorical_irls"
  )
  expect_identical(
    .scdesign3_scglm_method("auto", "poisson", ~ x, dat, scglm_fit = "approximate"),
    "newton_stein"
  )
  expect_identical(
    .scdesign3_scglm_method("auto", "poisson", ~ x, dat, scglm_fit = "exact"),
    "irls"
  )
  expect_identical(
    .scdesign3_scglm_method("auto", "nb", ~ x, dat, scglm_fit = "approximate"),
    "newton_stein"
  )
  expect_identical(
    .scdesign3_scglm_method("auto", "nb", ~ x, dat, scglm_fit = "exact"),
    "unsupported"
  )
})

test_that("fit_marginal auto avoids slow NB smooth scGLM path", {
  skip_if_not_installed("scGLM")
  skip_if_not_installed("mgcv")

  set.seed(102)
  n <- 60L
  g <- 3L
  dat <- data.frame(pseudotime = stats::runif(n))
  rownames(dat) <- paste0("cell", seq_len(n))
  count_mat <- matrix(stats::rnbinom(n * g, size = 3, mu = 3), nrow = n)
  rownames(count_mat) <- rownames(dat)
  colnames(count_mat) <- paste0("gene", seq_len(g))

  fit <- fit_marginal(
    data = list(count_mat = count_mat, dat = dat, filtered_gene = NULL),
    mu_formula = "s(pseudotime, bs = 'cr', k = 6)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    use_scglm = "auto"
  )

  expect_false(inherits(fit[[1]]$fit, "scdesign3_scglm"))
})

test_that("fit_marginal uses compressed categorical NBI backend for GAMLSS-style sigma models", {
  skip_if_not_installed("SingleCellExperiment")

  set.seed(103)
  n <- 120L
  g <- 4L
  group <- factor(rep(c("a", "b", "c"), length.out = n))
  dat <- data.frame(group = group)
  rownames(dat) <- paste0("cell", seq_len(n))
  group_id <- as.integer(group)
  mu <- matrix(c(2, 5, 9, 4, 7, 3, 6, 8, 10, 3, 5, 7), nrow = 3L)
  sigma <- matrix(c(0.2, 0.5, 0.8, 0.3, 0.7, 0.4, 0.5, 0.9, 0.6, 0.4, 0.8, 0.3), nrow = 3L)
  counts_cell_gene <- matrix(0, nrow = n, ncol = g)
  for (j in seq_len(g)) {
    counts_cell_gene[, j] <- stats::rnbinom(
      n,
      size = 1 / sigma[group_id, j],
      mu = mu[group_id, j]
    )
  }
  rownames(counts_cell_gene) <- rownames(dat)
  colnames(counts_cell_gene) <- paste0("gene", seq_len(g))
  counts_gene_cell <- t(counts_cell_gene)
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_gene_cell),
    colData = dat
  )

  fit <- fit_marginal(
    data = list(count_mat = counts_cell_gene, dat = dat, filtered_gene = NULL),
    mu_formula = "group",
    sigma_formula = "group",
    family_use = "nb",
    n_cores = 1,
    use_scglm = "auto"
  )

  expect_s3_class(fit[[1]]$fit, "scdesign3_categorical_nbi")
  expect_identical(fit[[1]]$fit$backend, "scDesign3::categorical_nbi_group_mle")
  expect_true(is.finite(stats::AIC(fit[[1]]$fit)))

  para <- extract_para(
    sce = sce,
    assay_use = "counts",
    marginal_list = fit,
    n_cores = 1,
    family_use = "nb",
    new_covariate = dat,
    data = dat
  )
  expected_means <- vapply(seq_len(g), function(j) {
    ave(counts_cell_gene[, j], group, FUN = mean)
  }, numeric(n))
  colnames(expected_means) <- colnames(counts_cell_gene)
  rownames(expected_means) <- rownames(counts_cell_gene)

  expect_equal(para$mean_mat, expected_means, tolerance = 1e-10)
  expect_true(all(is.finite(para$sigma_mat)))
  expect_true(all(para$sigma_mat > 0))
})

test_that("compressed categorical backend supports non-NB scDesign3 families", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("gamlss.dist")

  set.seed(104)
  n <- 90L
  g <- 3L
  group <- factor(rep(c("a", "b", "c"), length.out = n))
  group_id <- as.integer(group)
  dat <- data.frame(group = group)
  rownames(dat) <- paste0("cell", seq_len(n))

  make_sce <- function(counts_cell_gene) {
    rownames(counts_cell_gene) <- rownames(dat)
    colnames(counts_cell_gene) <- paste0("gene", seq_len(ncol(counts_cell_gene)))
    SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = t(counts_cell_gene)),
      colData = dat
    )
  }

  family_counts <- list(
    poisson = matrix(stats::rpois(n * g, lambda = rep(c(2, 5, 8), length.out = n)), nrow = n),
    binomial = matrix(stats::rbinom(n * g, size = 1, prob = rep(c(0.2, 0.5, 0.8), length.out = n)), nrow = n),
    gaussian = matrix(stats::rnorm(n * g, mean = rep(c(1, 4, 7), length.out = n), sd = rep(c(0.3, 0.8, 1.2), length.out = n)), nrow = n),
    zip = matrix(gamlss.dist::rZIP(n * g, mu = rep(c(2, 5, 8), length.out = n), sigma = rep(c(0.1, 0.3, 0.5), length.out = n)), nrow = n),
    zinb = matrix(gamlss.dist::rZINBI(n * g, mu = rep(c(2, 5, 8), length.out = n), sigma = rep(c(0.2, 0.5, 0.8), length.out = n), nu = rep(c(0.1, 0.2, 0.4), length.out = n)), nrow = n)
  )

  for (family_key in names(family_counts)) {
    counts <- family_counts[[family_key]]
    rownames(counts) <- rownames(dat)
    colnames(counts) <- paste0("gene", seq_len(ncol(counts)))
    sce <- make_sce(counts)
    fit <- fit_marginal(
      data = list(count_mat = counts, dat = dat, filtered_gene = NULL),
      mu_formula = "group",
      sigma_formula = "group",
      family_use = family_key,
      n_cores = 1,
      use_scglm = "auto"
    )

    expect_true(.scdesign3_is_categorical_dist(fit[[1]]$fit), info = family_key)
    expect_true(is.finite(stats::AIC(fit[[1]]$fit)), info = family_key)

    para <- extract_para(
      sce = sce,
      assay_use = "counts",
      marginal_list = fit,
      n_cores = 1,
      family_use = family_key,
      new_covariate = dat,
      data = dat
    )
    expect_true(all(is.finite(para$mean_mat)), info = family_key)
    if (family_key %in% c("gaussian", "nb", "zinb")) {
      expect_true(all(is.finite(para$sigma_mat)), info = family_key)
      expect_true(all(para$sigma_mat > 0), info = family_key)
    }
    if (family_key %in% c("zip", "zinb")) {
      expect_true(all(is.finite(as.matrix(para$zero_mat))), info = family_key)
      expect_true(all(as.matrix(para$zero_mat) >= 0 & as.matrix(para$zero_mat) <= 1), info = family_key)
    }

    u <- convert_u(
      sce = sce,
      assay_use = "counts",
      marginal_list = fit,
      data = dat,
      DT = !identical(family_key, "gaussian"),
      pseudo_obs = FALSE,
      epsilon = 1e-6,
      n_cores = 1,
      family_use = family_key,
      parallelization = "mcmapply",
      BPPARAM = NULL
    )
    expect_equal(dim(u), c(n, g), info = family_key)
    expect_true(all(is.finite(u)), info = family_key)
    expect_true(all(u > 0 & u < 1), info = family_key)
  }
})

test_that("perform_lrt accepts categorical fast-path models", {
  set.seed(105)
  n <- 90L
  g <- 4L
  group <- factor(rep(c("a", "b", "c"), length.out = n))
  dat <- data.frame(group = group)
  rownames(dat) <- paste0("cell", seq_len(n))
  counts <- matrix(
    stats::rpois(n * g, lambda = rep(c(2, 5, 8), length.out = n)),
    nrow = n
  )
  rownames(counts) <- rownames(dat)
  colnames(counts) <- paste0("gene", seq_len(g))

  alt <- fit_marginal(
    data = list(count_mat = counts, dat = dat, filtered_gene = NULL),
    mu_formula = "group",
    sigma_formula = "group",
    family_use = "poisson",
    n_cores = 1,
    use_scglm = "auto"
  )
  null <- fit_marginal(
    data = list(count_mat = counts, dat = dat, filtered_gene = NULL),
    mu_formula = "1",
    sigma_formula = "1",
    family_use = "poisson",
    n_cores = 1,
    use_scglm = "auto"
  )

  lrt <- perform_lrt(
    lapply(alt, function(x) x$fit),
    lapply(null, function(x) x$fit)
  )

  expect_equal(nrow(lrt), g)
  expect_true(all(is.finite(as.numeric(lrt$p_value)) | is.na(lrt$p_value)))
})

test_that("shared fixed-smooth backend matches original mgcv path", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("mgcv")

  set.seed(106)
  n <- 90L
  g <- 3L
  dat <- data.frame(x = stats::runif(n))
  rownames(dat) <- paste0("cell", seq_len(n))
  eta <- 0.3 + sin(2 * pi * dat$x)
  forms <- c(
    "s(x, bs = 'tp', k = 8, fx = TRUE)",
    "s(x, bs = 'tp', k = 8, sp = 0.4)"
  )
  families <- c("poisson", "gaussian", "binomial")

  for (family_key in families) {
    counts <- switch(
      family_key,
      poisson = matrix(stats::rpois(n * g, lambda = rep(exp(eta), g)), nrow = n),
      gaussian = matrix(stats::rnorm(n * g, mean = rep(eta, g), sd = 0.4), nrow = n),
      binomial = matrix(stats::rbinom(n * g, size = 1, prob = stats::plogis(rep(eta, g))), nrow = n)
    )
    rownames(counts) <- rownames(dat)
    colnames(counts) <- paste0("gene", seq_len(g))
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = t(counts)),
      colData = dat
    )

    for (mu_formula in forms) {
      fast <- fit_marginal(
        data = list(count_mat = counts, dat = dat, filtered_gene = NULL),
        mu_formula = mu_formula,
        sigma_formula = "1",
        family_use = family_key,
        n_cores = 2,
        use_scglm = "auto",
        scglm_batch_size = 1L
      )
      slow <- fit_marginal(
        data = list(count_mat = counts, dat = dat, filtered_gene = NULL),
        mu_formula = mu_formula,
        sigma_formula = "1",
        family_use = family_key,
        n_cores = 1,
        use_scglm = "never"
      )

      expect_s3_class(fast[[1]]$fit, "scdesign3_scglm")
      expect_s3_class(slow[[1]]$fit, "gam")

      fast_para <- extract_para(
        sce = sce,
        assay_use = "counts",
        marginal_list = fast,
        n_cores = 1,
        family_use = family_key,
        new_covariate = dat,
        data = dat
      )
      slow_para <- extract_para(
        sce = sce,
        assay_use = "counts",
        marginal_list = slow,
        n_cores = 1,
        family_use = family_key,
        new_covariate = dat,
        data = dat
      )

      expect_equal(fast_para$mean_mat, slow_para$mean_mat, tolerance = 1e-8, info = paste(family_key, mu_formula))
      if (identical(family_key, "gaussian")) {
        expect_equal(fast_para$sigma_mat, slow_para$sigma_mat, tolerance = 1e-8, info = paste(family_key, mu_formula))
      }
      expect_equal(
        stats::AIC(fast[[1]]$fit),
        stats::AIC(slow[[1]]$fit),
        tolerance = 1e-8,
        info = paste(family_key, mu_formula)
      )
    }
  }
})

test_that("automatic penalized smooths fall back to original mgcv path", {
  skip_if_not_installed("mgcv")

  set.seed(107)
  n <- 60L
  g <- 2L
  dat <- data.frame(x = stats::runif(n))
  rownames(dat) <- paste0("cell", seq_len(n))
  mu <- exp(0.2 + sin(2 * pi * dat$x))
  counts <- matrix(stats::rpois(n * g, lambda = rep(mu, g)), nrow = n)
  rownames(counts) <- rownames(dat)
  colnames(counts) <- paste0("gene", seq_len(g))

  fit <- fit_marginal(
    data = list(count_mat = counts, dat = dat, filtered_gene = NULL),
    mu_formula = "s(x, bs = 'tp', k = 8)",
    sigma_formula = "1",
    family_use = "poisson",
    n_cores = 1,
    use_scglm = "auto"
  )

  expect_s3_class(fit[[1]]$fit, "gam")
  expect_false(inherits(fit[[1]]$fit, "scdesign3_scglm"))
})
