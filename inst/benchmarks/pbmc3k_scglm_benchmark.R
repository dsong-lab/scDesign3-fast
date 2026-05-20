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
n_genes <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1500L
n_clusters <- if (length(args) >= 3L) as.integer(args[[3L]]) else 8L
out_rds <- if (length(args) >= 4L) args[[4L]] else "/tmp/pbmc3k_scglm_benchmark.rds"

read_tsv_first <- function(path) {
  utils::read.delim(path, header = FALSE, stringsAsFactors = FALSE)
}

counts <- Matrix::readMM(file.path(pbmc_dir, "matrix.mtx"))
genes <- read_tsv_first(file.path(pbmc_dir, "genes.tsv"))
barcodes <- read_tsv_first(file.path(pbmc_dir, "barcodes.tsv"))[[1L]]
gene_names <- make.unique(genes[[2L]])
rownames(counts) <- gene_names
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
set.seed(1)
cluster <- factor(stats::kmeans(pc$x, centers = n_clusters, iter.max = 50L)$cluster)

count_mat <- t(as.matrix(counts_sel))
dat <- data.frame(cluster = cluster)
rownames(dat) <- rownames(count_mat)
bench_data <- list(count_mat = count_mat, dat = dat, filtered_gene = NULL)

gc()
scglm_time <- system.time({
  fit_scglm <- fit_marginal(
    data = bench_data,
    mu_formula = "cluster",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    use_scglm = "auto"
  )
})

gc()
mgcv_time <- system.time({
  fit_mgcv <- fit_marginal(
    data = bench_data,
    mu_formula = "cluster",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    use_scglm = "never"
  )
})

pred_scglm <- unlist(lapply(fit_scglm, function(z) as.numeric(stats::predict(z$fit, type = "response"))), use.names = FALSE)
pred_mgcv <- unlist(lapply(fit_mgcv, function(z) as.numeric(stats::predict(z$fit, type = "response"))), use.names = FALSE)
theta_scglm <- vapply(fit_scglm, function(z) z$fit$theta, numeric(1))
theta_mgcv <- vapply(fit_mgcv, function(z) z$fit$family$getTheta(TRUE), numeric(1))
aic_scglm <- vapply(fit_scglm, function(z) stats::AIC(z$fit), numeric(1))
aic_mgcv <- vapply(fit_mgcv, function(z) stats::AIC(z$fit), numeric(1))

summary <- data.frame(
  dataset = "PBMC3k 10x",
  cells = nrow(count_mat),
  genes = ncol(count_mat),
  clusters = nlevels(cluster),
  formula = "gene ~ cluster",
  family = "nb",
  scglm_backend = fit_scglm[[1L]]$fit$backend,
  scglm_elapsed = unname(scglm_time[["elapsed"]]),
  mgcv_elapsed = unname(mgcv_time[["elapsed"]]),
  speedup = unname(mgcv_time[["elapsed"]] / scglm_time[["elapsed"]]),
  pred_cor = stats::cor(pred_scglm, pred_mgcv),
  pred_rmse = sqrt(mean((pred_scglm - pred_mgcv)^2)),
  theta_cor = stats::cor(theta_scglm, theta_mgcv),
  theta_rmse = sqrt(mean((theta_scglm - theta_mgcv)^2)),
  median_abs_aic_diff = stats::median(abs(aic_scglm - aic_mgcv)),
  max_abs_aic_diff = max(abs(aic_scglm - aic_mgcv)),
  scglm_object_mb = as.numeric(utils::object.size(fit_scglm)) / 1024^2,
  mgcv_object_mb = as.numeric(utils::object.size(fit_mgcv)) / 1024^2,
  stringsAsFactors = FALSE
)

saveRDS(
  list(
    summary = summary,
    theta = data.frame(gene = colnames(count_mat), scglm = theta_scglm, mgcv = theta_mgcv),
    aic = data.frame(gene = colnames(count_mat), scglm = aic_scglm, mgcv = aic_mgcv)
  ),
  out_rds
)

print(summary)
cat("Saved benchmark result:", out_rds, "\n")
