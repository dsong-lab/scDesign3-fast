#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L) {
  stop("Usage: scdesign3_celltype_vignette_original_vs_new.R <scdesign3_lib> <extra_lib> <label> <sce_rds> <out_prefix>")
}

scdesign3_lib <- args[[1L]]
extra_lib <- args[[2L]]
label <- args[[3L]]
sce_rds <- args[[4L]]
out_prefix <- args[[5L]]

.libPaths(c(scdesign3_lib, extra_lib, .libPaths()))

suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

rss_mb <- function() {
  out <- suppressWarnings(system(sprintf("ps -o rss= -p %s", Sys.getpid()), intern = TRUE))
  if (length(out) == 0L) return(NA_real_)
  as.numeric(trimws(out[[1L]])) / 1024
}

obj_mb <- function(x) {
  as.numeric(object.size(x)) / 1024^2
}

time_step <- function(name, fun) {
  invisible(gc())
  rss_start <- rss_mb()
  t0 <- proc.time()[["elapsed"]]
  value <- fun()
  elapsed <- proc.time()[["elapsed"]] - t0
  invisible(gc())
  list(
    value = value,
    timing = data.frame(
      label = label,
      step = name,
      elapsed_sec = elapsed,
      rss_mb_after = rss_mb(),
      rss_mb_delta = rss_mb() - rss_start,
      stringsAsFactors = FALSE
    )
  )
}

make_modified_celltype_covariate <- function(input_data, n_cell) {
  ct_levels <- levels(as.factor(input_data$cell_type))
  ct_prop <- c(0, 0, 0.2, 0.8)
  if (length(ct_levels) != length(ct_prop)) {
    stop("Unexpected number of cell types in vignette dataset.")
  }
  n_per_type <- as.integer(round(ct_prop * n_cell))
  n_per_type[length(n_per_type)] <- n_cell - sum(n_per_type[-length(n_per_type)])
  cell_type <- rep(ct_levels, n_per_type)
  data.frame(
    cell_type = factor(cell_type, levels = ct_levels),
    corr_group = cell_type,
    row.names = paste0("Cell", seq_along(cell_type))
  )
}

set.seed(123)
sce <- readRDS(sce_rds)

family_use <- "nb"
n_cores_fit <- 2L
mu_formula <- "cell_type"
sigma_formula <- "cell_type"

construct_res <- time_step("construct_data", function() construct_data(
  sce = sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = NULL,
  corr_by = "cell_type"
))
constructed <- construct_res$value

marginal_res <- time_step("fit_marginal", function() fit_marginal(
  data = constructed,
  predictor = "gene",
  mu_formula = mu_formula,
  sigma_formula = sigma_formula,
  family_use = family_use,
  n_cores = n_cores_fit,
  usebam = FALSE,
  parallelization = "pbmcmapply"
))
marginal_list <- marginal_res$value

set.seed(123)
copula_res <- time_step("fit_copula", function() fit_copula(
  sce = sce,
  assay_use = "counts",
  input_data = constructed$dat,
  marginal_list = marginal_list,
  family_use = family_use,
  copula = "gaussian",
  n_cores = n_cores_fit
))
copula_list <- copula_res$value

extract_same_res <- time_step("extract_para_same", function() extract_para(
  sce = sce,
  assay_use = "counts",
  marginal_list = marginal_list,
  n_cores = 1L,
  family_use = family_use,
  new_covariate = NULL,
  data = constructed$dat
))
para_same <- extract_same_res$value

set.seed(123)
simu_same_res <- time_step("simu_new_same", function() simu_new(
  sce = sce,
  assay_use = "counts",
  mean_mat = para_same$mean_mat,
  sigma_mat = para_same$sigma_mat,
  zero_mat = para_same$zero_mat,
  quantile_mat = NULL,
  copula_list = copula_list$copula_list,
  n_cores = 1L,
  family_use = family_use,
  input_data = constructed$dat,
  new_covariate = constructed$newCovariate,
  important_feature = copula_list$important_feature,
  filtered_gene = constructed$filtered_gene
))
new_count_same <- simu_same_res$value

new_ct <- make_modified_celltype_covariate(constructed$dat, ncol(sce))

extract_modified_res <- time_step("extract_para_modified", function() extract_para(
  sce = sce,
  assay_use = "counts",
  marginal_list = marginal_list,
  n_cores = 1L,
  family_use = family_use,
  new_covariate = new_ct,
  data = constructed$dat
))
para_modified <- extract_modified_res$value

set.seed(123)
simu_modified_res <- time_step("simu_new_modified", function() simu_new(
  sce = sce,
  assay_use = "counts",
  mean_mat = para_modified$mean_mat,
  sigma_mat = para_modified$sigma_mat,
  zero_mat = para_modified$zero_mat,
  quantile_mat = NULL,
  copula_list = copula_list$copula_list,
  n_cores = 1L,
  family_use = family_use,
  input_data = constructed$dat,
  new_covariate = new_ct,
  important_feature = copula_list$important_feature,
  filtered_gene = constructed$filtered_gene
))
new_count_modified <- simu_modified_res$value

timing <- do.call(
  rbind,
  list(
    construct_res$timing,
    marginal_res$timing,
    copula_res$timing,
    extract_same_res$timing,
    simu_same_res$timing,
    extract_modified_res$timing,
    simu_modified_res$timing
  )
)
timing <- rbind(
  timing,
  data.frame(
    label = label,
    step = "total",
    elapsed_sec = sum(timing$elapsed_sec),
    rss_mb_after = rss_mb(),
    rss_mb_delta = NA_real_,
    stringsAsFactors = FALSE
  )
)

count_mat <- constructed$count_mat
same_dense <- as.matrix(new_count_same)
modified_dense <- as.matrix(new_count_modified)
summary <- data.frame(
  label = label,
  package_version = as.character(utils::packageVersion("scDesign3")),
  n_genes = nrow(sce),
  n_cells = ncol(sce),
  filtered_genes = length(constructed$filtered_gene),
  count_mat_class = paste(class(count_mat), collapse = "/"),
  count_mat_mb = obj_mb(count_mat),
  marginal_list_mb = obj_mb(marginal_list),
  copula_list_mb = obj_mb(copula_list),
  para_same_mb = obj_mb(para_same),
  para_modified_mb = obj_mb(para_modified),
  new_count_same_mb = obj_mb(new_count_same),
  new_count_modified_mb = obj_mb(new_count_modified),
  input_zero_prop = mean(SummarizedExperiment::assay(sce, "counts") == 0),
  same_zero_prop = mean(same_dense == 0),
  modified_zero_prop = mean(modified_dense == 0),
  same_mean = mean(same_dense),
  modified_mean = mean(modified_dense),
  same_sum = sum(same_dense),
  modified_sum = sum(modified_dense),
  stringsAsFactors = FALSE
)

out <- list(
  label = label,
  timing = timing,
  summary = summary,
  filtered_gene = constructed$filtered_gene,
  new_ct = new_ct,
  mean_mat_same = para_same$mean_mat,
  sigma_mat_same = para_same$sigma_mat,
  zero_mat_same = para_same$zero_mat,
  mean_mat_modified = para_modified$mean_mat,
  sigma_mat_modified = para_modified$sigma_mat,
  zero_mat_modified = para_modified$zero_mat,
  new_count_same = new_count_same,
  new_count_modified = new_count_modified
)

saveRDS(out, paste0(out_prefix, ".rds"))
utils::write.csv(timing, paste0(out_prefix, "_timing.csv"), row.names = FALSE)
utils::write.csv(summary, paste0(out_prefix, "_summary.csv"), row.names = FALSE)

print(timing)
print(summary)
