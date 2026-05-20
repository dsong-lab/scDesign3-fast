#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4L) {
  stop("Usage: scdesign3_vignette_original_vs_new.R <lib_path> <label> <sce_rds> <out_prefix>")
}

lib_path <- args[[1L]]
label <- args[[2L]]
sce_rds <- args[[3L]]
out_prefix <- args[[4L]]

.libPaths(c(lib_path, .libPaths()))

suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
})

rss_mb <- function() {
  out <- suppressWarnings(system(sprintf("ps -o rss= -p %s", Sys.getpid()), intern = TRUE))
  if (length(out) == 0L) {
    return(NA_real_)
  }
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

set.seed(123)
sce_full <- readRDS(sce_rds)
sce <- sce_full[seq_len(100L), ]

family_use <- "nb"
n_cores <- 2L
mu_formula <- "s(pseudotime, k = 10, bs = 'cr')"
sigma_formula <- "s(pseudotime, k = 5, bs = 'cr')"
important_feature <- "all"

construct_res <- time_step("construct_data", function() construct_data(
  sce = sce,
  assay_use = "counts",
  celltype = "cell_type",
  pseudotime = "pseudotime",
  spatial = NULL,
  other_covariates = NULL,
  corr_by = "1"
))
constructed <- construct_res$value

marginal_res <- time_step("fit_marginal", function() fit_marginal(
  data = constructed,
  mu_formula = mu_formula,
  sigma_formula = sigma_formula,
  family_use = family_use,
  n_cores = n_cores,
  parallelization = "mcmapply"
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
  empirical_quantile = FALSE,
  if_sparse = TRUE,
  important_feature = important_feature,
  n_cores = n_cores,
  parallelization = "mcmapply"
))
copula_list <- copula_res$value

extract_res <- time_step("extract_para", function() extract_para(
  sce = sce,
  assay_use = "counts",
  marginal_list = marginal_list,
  n_cores = 1L,
  family_use = family_use,
  new_covariate = constructed$newCovariate,
  data = constructed$dat,
  parallelization = "mcmapply"
))
para <- extract_res$value

set.seed(123)
simu_res <- time_step("simu_new", function() simu_new(
  sce = sce,
  assay_use = "counts",
  mean_mat = para$mean_mat,
  sigma_mat = para$sigma_mat,
  zero_mat = para$zero_mat,
  quantile_mat = NULL,
  copula_list = copula_list$copula_list,
  n_cores = 1L,
  fastmvn = FALSE,
  family_use = family_use,
  input_data = constructed$dat,
  new_covariate = constructed$newCovariate,
  important_feature = copula_list$important_feature,
  parallelization = "mcmapply",
  filtered_gene = constructed$filtered_gene
))
new_count <- simu_res$value

timing <- do.call(
  rbind,
  list(
    construct_res$timing,
    marginal_res$timing,
    copula_res$timing,
    extract_res$timing,
    simu_res$timing
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
new_count_dense <- as.matrix(new_count)
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
  para_mb = obj_mb(para),
  new_count_mb = obj_mb(new_count),
  input_zero_prop = mean(SummarizedExperiment::assay(sce, "counts") == 0),
  simulated_zero_prop = mean(new_count_dense == 0),
  simulated_mean = mean(new_count_dense),
  simulated_sum = sum(new_count_dense),
  stringsAsFactors = FALSE
)

out <- list(
  label = label,
  timing = timing,
  summary = summary,
  filtered_gene = constructed$filtered_gene,
  mean_mat = para$mean_mat,
  sigma_mat = para$sigma_mat,
  zero_mat = para$zero_mat,
  new_count = new_count
)

saveRDS(out, paste0(out_prefix, ".rds"))
utils::write.csv(timing, paste0(out_prefix, "_timing.csv"), row.names = FALSE)
utils::write.csv(summary, paste0(out_prefix, "_summary.csv"), row.names = FALSE)

print(timing)
print(summary)
