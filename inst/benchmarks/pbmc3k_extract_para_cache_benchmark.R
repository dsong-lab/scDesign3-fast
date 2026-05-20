#!/usr/bin/env Rscript

bench_libs <- Sys.getenv("SCGLM_BENCH_LIBS", "")
if (nzchar(bench_libs)) {
  .libPaths(c(strsplit(bench_libs, .Platform$path.sep, fixed = TRUE)[[1L]], .libPaths()))
}

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(SingleCellExperiment)
  library(scDesign3)
})

args <- commandArgs(trailingOnly = TRUE)
pbmc_dir <- if (length(args) >= 1L) args[[1L]] else "/home/dosong/scalableGLM/pbmc3k_data/filtered_gene_bc_matrices/hg19"
n_genes <- if (length(args) >= 2L) as.integer(args[[2L]]) else 500L
k <- if (length(args) >= 3L) as.integer(args[[3L]]) else 10L
n_new <- if (length(args) >= 4L) as.integer(args[[4L]]) else 2700L
out_rds <- if (length(args) >= 5L) args[[5L]] else "/tmp/pbmc3k_extract_para_cache_benchmark.rds"

read_tsv_first <- function(path) {
  utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE)
}

counts <- Matrix::readMM(file.path(pbmc_dir, "matrix.mtx"))
genes <- read_tsv_first(file.path(pbmc_dir, "genes.tsv"))
barcodes <- read_tsv_first(file.path(pbmc_dir, "barcodes.tsv"))[[1L]]
rownames(counts) <- make.unique(genes[[2L]])
colnames(counts) <- barcodes

detected <- Matrix::rowSums(counts > 0)
total <- Matrix::rowSums(counts)
keep <- detected >= max(10L, floor(0.01 * ncol(counts))) & total > 0
counts <- counts[keep, , drop = FALSE]
gene_score <- Matrix::rowMeans(counts) * log1p(detected[keep])
sel <- order(gene_score, decreasing = TRUE)[seq_len(min(n_genes, length(gene_score)))]
counts_sel <- counts[sel, , drop = FALSE]

norm <- t(t(counts_sel) / pmax(Matrix::colSums(counts_sel), 1)) * 1e4
log_norm <- log1p(norm)
gene_var <- matrixStats::rowVars(as.matrix(log_norm))
pca_genes <- order(gene_var, decreasing = TRUE)[seq_len(min(1000L, length(gene_var)))]
x_pca <- t(as.matrix(log_norm[pca_genes, , drop = FALSE]))
x_pca <- scale(x_pca, center = TRUE, scale = FALSE)
pc <- irlba::prcomp_irlba(x_pca, n = min(20L, ncol(x_pca) - 1L), center = FALSE, scale. = FALSE)
pseudotime <- as.numeric(scale(pc$x[, 1L]))

count_mat <- t(as.matrix(counts_sel))
dat <- data.frame(pseudotime = pseudotime)
rownames(dat) <- rownames(count_mat)
bench_data <- list(count_mat = count_mat, dat = dat, filtered_gene = NULL)
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = t(count_mat))
)
mu_formula <- paste0("s(pseudotime, bs = 'cr', k = ", k, ")")

gc()
fit_time <- system.time({
  fit <- fit_marginal(
    data = bench_data,
    mu_formula = mu_formula,
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    use_scglm = "never"
  )
})

new_covariate <- data.frame(
  pseudotime = seq(min(dat$pseudotime), max(dat$pseudotime), length.out = n_new)
)
rownames(new_covariate) <- paste0("newcell", seq_len(nrow(new_covariate)))

gc()
cached_extract_time <- system.time({
  para_cached <- extract_para(
    sce = sce,
    marginal_list = fit,
    n_cores = 1,
    family_use = "nb",
    new_covariate = new_covariate,
    data = dat
  )
})

gc()
direct_predict_time <- system.time({
  direct_mean <- sapply(fit, function(x) {
    stats::predict(x$fit, type = "response", newdata = new_covariate)
  })
  direct_sigma <- sapply(fit, function(x) {
    rep(1 / x$fit$family$getTheta(TRUE), nrow(new_covariate))
  })
})

summary <- data.frame(
  dataset = "PBMC3k 10x",
  cells_fit = nrow(count_mat),
  cells_predict = nrow(new_covariate),
  genes = ncol(count_mat),
  formula = paste0("gene ~ ", mu_formula),
  family = "nb",
  fit_elapsed = unname(fit_time[["elapsed"]]),
  cached_extract_elapsed = unname(cached_extract_time[["elapsed"]]),
  direct_predict_elapsed = unname(direct_predict_time[["elapsed"]]),
  extract_speedup = unname(direct_predict_time[["elapsed"]] / cached_extract_time[["elapsed"]]),
  pred_cor = stats::cor(as.vector(para_cached$mean_mat), as.vector(direct_mean)),
  pred_rmse = sqrt(mean((as.vector(para_cached$mean_mat) - as.vector(direct_mean))^2)),
  sigma_rmse = sqrt(mean((as.vector(para_cached$sigma_mat) - as.vector(direct_sigma))^2)),
  cached_object_mb = as.numeric(utils::object.size(para_cached)) / 1024^2,
  direct_object_mb = as.numeric(utils::object.size(list(mean = direct_mean, sigma = direct_sigma))) / 1024^2,
  stringsAsFactors = FALSE
)

saveRDS(
  list(
    summary = summary,
    mean_diff_summary = summary(as.vector(para_cached$mean_mat) - as.vector(direct_mean))
  ),
  out_rds
)

print(summary)
cat("Saved benchmark result:", out_rds, "\n")
