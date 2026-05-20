test_that("construct_data preserves sparse count matrices", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Matrix")

  counts <- matrix(
    c(0, 1, 0, 0, 2, 0, 0, 0, 3, 0, 0, 4),
    nrow = 3,
    byrow = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(counts, sparse = TRUE)),
    colData = data.frame(group = factor(c("a", "a", "b", "b")), row.names = colnames(counts))
  )

  out <- construct_data(
    sce = sce,
    assay_use = "counts",
    celltype = "group",
    pseudotime = NULL,
    spatial = NULL,
    other_covariates = NULL,
    corr_by = "1"
  )

  expect_s4_class(out$count_mat, "sparseMatrix")
  expect_equal(as.matrix(out$count_mat), t(counts))
})

test_that("simu_new uses aligned blockwise quantiles without full gene expansion", {
  skip_if_not_installed("SingleCellExperiment")

  n <- 5L
  g <- 4L
  gene_names <- paste0("gene", seq_len(g))
  cell_names <- paste0("cell", seq_len(n))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = g, ncol = n, dimnames = list(gene_names, cell_names)))
  )
  mean_mat <- matrix(2, nrow = n, ncol = g, dimnames = list(cell_names, gene_names))
  mean_mat[2, "gene2"] <- 0
  sigma_mat <- matrix(1, nrow = n, ncol = g, dimnames = dimnames(mean_mat))
  zero_mat <- matrix(0, nrow = n, ncol = g, dimnames = dimnames(mean_mat))
  quantile_mat <- matrix(
    seq(0.1, 0.9, length.out = n * 2L),
    nrow = n,
    dimnames = list(cell_names, c("gene1", "gene2"))
  )
  input_data <- data.frame(corr_group = rep(1, n), row.names = cell_names)

  out <- simu_new(
    sce = sce,
    mean_mat = mean_mat,
    sigma_mat = sigma_mat,
    zero_mat = zero_mat,
    quantile_mat = quantile_mat,
    copula_list = NULL,
    n_cores = 1,
    family_use = "poisson",
    input_data = input_data,
    new_covariate = input_data,
    filtered_gene = c("gene3", "gene4")
  )

  expected <- matrix(0, nrow = g, ncol = n, dimnames = list(gene_names, cell_names))
  expected["gene1", ] <- stats::qpois(quantile_mat[, "gene1"], lambda = mean_mat[, "gene1"])
  expected["gene2", ] <- stats::qpois(quantile_mat[, "gene2"], lambda = mean_mat[, "gene2"])
  expected["gene2", mean_mat[, "gene2"] == 0] <- 0

  expect_equal(out, expected)
})

test_that("simu_new falls back to positional quantile alignment when names do not match", {
  skip_if_not_installed("SingleCellExperiment")

  n <- 4L
  g <- 3L
  gene_names <- paste0("gene", seq_len(g))
  cell_names <- paste0("cell", seq_len(n))
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = g, ncol = n, dimnames = list(gene_names, cell_names)))
  )
  mean_mat <- matrix(3, nrow = n, ncol = g, dimnames = list(cell_names, gene_names))
  sigma_mat <- matrix(1, nrow = n, ncol = g, dimnames = dimnames(mean_mat))
  zero_mat <- matrix(0, nrow = n, ncol = g, dimnames = dimnames(mean_mat))
  quantile_mat <- matrix(
    c(0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9),
    nrow = n,
    dimnames = list(as.character(seq_len(n)), c("x", "y"))
  )
  input_data <- data.frame(corr_group = rep(1, n), row.names = cell_names)

  out <- simu_new(
    sce = sce,
    mean_mat = mean_mat,
    sigma_mat = sigma_mat,
    zero_mat = zero_mat,
    quantile_mat = quantile_mat,
    copula_list = NULL,
    n_cores = 1,
    family_use = "poisson",
    input_data = input_data,
    new_covariate = input_data,
    filtered_gene = "gene3"
  )

  expected <- matrix(0, nrow = g, ncol = n, dimnames = list(gene_names, cell_names))
  expected["gene1", ] <- stats::qpois(quantile_mat[, 1], lambda = mean_mat[, 1])
  expected["gene2", ] <- stats::qpois(quantile_mat[, 2], lambda = mean_mat[, 2])

  expect_equal(out, expected)
})
