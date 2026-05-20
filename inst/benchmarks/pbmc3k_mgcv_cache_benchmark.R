#!/usr/bin/env Rscript

bench_libs <- Sys.getenv("SCGLM_BENCH_LIBS", "")
if (nzchar(bench_libs)) {
  .libPaths(c(strsplit(bench_libs, .Platform$path.sep, fixed = TRUE)[[1L]], .libPaths()))
}

suppressPackageStartupMessages({
  library(Matrix)
  library(irlba)
  library(scDesign3)
})

args <- commandArgs(trailingOnly = TRUE)
pbmc_dir <- if (length(args) >= 1L) args[[1L]] else "/home/dosong/scalableGLM/pbmc3k_data/filtered_gene_bc_matrices/hg19"
n_genes <- if (length(args) >= 2L) as.integer(args[[2L]]) else 500L
k <- if (length(args) >= 3L) as.integer(args[[3L]]) else 10L
out_rds <- if (length(args) >= 4L) args[[4L]] else "/tmp/pbmc3k_mgcv_cache_benchmark.rds"

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
mu_formula <- paste0("s(pseudotime, bs = 'cr', k = ", k, ")")
direct_formula <- stats::as.formula(paste0("gene ~ ", mu_formula))

gc()
cached_time <- system.time({
  fit_cached <- fit_marginal(
    data = bench_data,
    mu_formula = mu_formula,
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    use_scglm = "never"
  )
})

gc()
direct_time <- system.time({
  fit_direct <- lapply(seq_len(ncol(count_mat)), function(j) {
    dat_j <- dat
    dat_j$gene <- count_mat[, j]
    mgcv::gam(direct_formula, data = dat_j, family = "nb", method = "REML")
  })
})

pred_cached <- unlist(lapply(fit_cached, function(z) as.numeric(stats::predict(z$fit, type = "response"))), use.names = FALSE)
pred_direct <- unlist(lapply(fit_direct, function(z) as.numeric(stats::predict(z, type = "response"))), use.names = FALSE)
theta_cached <- vapply(fit_cached, function(z) z$fit$family$getTheta(TRUE), numeric(1))
theta_direct <- vapply(fit_direct, function(z) z$family$getTheta(TRUE), numeric(1))
aic_cached <- vapply(fit_cached, function(z) stats::AIC(z$fit), numeric(1))
aic_direct <- vapply(fit_direct, function(z) stats::AIC(z), numeric(1))

summary <- data.frame(
  dataset = "PBMC3k 10x",
  cells = nrow(count_mat),
  genes = ncol(count_mat),
  formula = paste0("gene ~ ", mu_formula),
  family = "nb",
  cached_elapsed = unname(cached_time[["elapsed"]]),
  direct_mgcv_elapsed = unname(direct_time[["elapsed"]]),
  speedup = unname(direct_time[["elapsed"]] / cached_time[["elapsed"]]),
  pred_cor = stats::cor(pred_cached, pred_direct),
  pred_rmse = sqrt(mean((pred_cached - pred_direct)^2)),
  theta_cor = stats::cor(theta_cached, theta_direct),
  theta_rmse = sqrt(mean((theta_cached - theta_direct)^2)),
  median_abs_aic_diff = stats::median(abs(aic_cached - aic_direct)),
  max_abs_aic_diff = max(abs(aic_cached - aic_direct)),
  cached_object_mb = as.numeric(utils::object.size(fit_cached)) / 1024^2,
  direct_object_mb = as.numeric(utils::object.size(fit_direct)) / 1024^2,
  stringsAsFactors = FALSE
)

saveRDS(
  list(
    summary = summary,
    theta = data.frame(gene = colnames(count_mat), cached = theta_cached, direct = theta_direct),
    aic = data.frame(gene = colnames(count_mat), cached = aic_cached, direct = aic_direct)
  ),
  out_rds
)

print(summary)
cat("Saved benchmark result:", out_rds, "\n")
