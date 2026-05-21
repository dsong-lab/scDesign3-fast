#' The wrapper for the whole scDesign3 pipeline
#'
#' \code{scdesign3} takes the input data, fits the model and
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assay_use A string which indicates the assay you will use in the sce. Default is 'counts'.
#' @param celltype A string of the name of cell type variable in the \code{colData} of the sce. Default is 'cell_type'.
#' @param pseudotime A string or a string vector of the name of pseudotime and (if exist)
#' multiple lineages. Default is NULL.
#' @param spatial A length two string vector of the names of spatial coordinates. Default is NULL.
#' @param other_covariates A string or a string vector of the other covariates you want to include in the data.
#' @param ncell The number of cell you want to simulate. Default is \code{dim(sce)[2]} (the same number as the input data).
#' @param mu_formula A string of the mu parameter formula
#' @param sigma_formula A string of the sigma parameter formula
#' @param family_use A string of the marginal distribution.
#' Must be one of 'poisson', 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param n_cores An integer. The number of cores to use.
#' @param correlation_function A string. If 'default', use the package default correlation implementation; if 'coop', use the implementation from \code{coop}, which calls BLAS.
#' @param usebam A logic variable. If use \code{\link[mgcv]{bam}} for acceleration in marginal fitting.
#' @param use_scglm A string indicating whether to use the optional batched
#' \code{scGLM} marginal backend. Must be one of \code{"auto"},
#' \code{"always"}, or \code{"never"}. For categorical shared designs,
#' \code{"auto"} uses the scGLM categorical closed-form backend by default.
#' @param scglm_method A string selecting the scGLM matrix backend.
#' @param scglm_batch_size Number of features to process per scGLM batch.
#' @param edf_flexible A logic variable. It is used for accelerating for spatial model if k is large in 'mu_formula'. Default is FALSE.
#' @param corr_formula A string of the correlation structure.
#' @param empirical_quantile Please only use it if you clearly know what will happen! A logic variable. If TRUE, DO NOT fit the copula and use the EMPIRICAL CDF values of the original data; it will make the simulated data fixed (no randomness). Default is FALSE. Only works if ncell is the same as your original data.
#' @param copula A string of the copula choice. Must be one of 'gaussian' or 'vine'. Default is 'gaussian'. Note that vine copula may have better modeling of high-dimensions, but can be very slow when features are >1000.
#' @param if_sparse A logic variable. Only works for Gaussian copula (\code{family_set = "gaussian"}). If TRUE, a thresholding strategy will make the corr matrix sparse.
#' @param fastmvn An logical variable. If TRUE, the sampling of multivariate Gaussian is done by \code{mvnfast}, otherwise by \code{mvtnorm}. Default is FALSE. It only matters for Gaussian copula.
#' @param DT A logic variable. If TRUE, perform the distributional transformation
#' to make the discrete data 'continuous'. This is useful for discrete distributions (e.g., Poisson, NB).
#' Default is TRUE. Note that for continuous data (e.g., Gaussian), DT does not make sense and should be set as FALSE.
#' @param pseudo_obs A logic variable. If TRUE, use the empirical quantiles instead of theoretical quantiles for fitting copula.
#' Default is FALSE.
#' @param family_set A string or a string vector of the bivariate copula families. Default is c("gauss", "indep"). For more information please check package \code{rvinecoplib}.
#' @param important_feature A numeric value or vector which indicates whether a gene will be used in correlation estimation or not. If this is a numeric value, then
#' gene with zero proportion greater than this value will be excluded form gene-gene correlation estimation. If this is a vector, then this should
#' be a logical vector with length equal to the number of genes in \code{sce}. \code{TRUE} in the logical vector means the corresponding gene will be included in
#' gene-gene correlation estimation and \code{FALSE} in the logical vector means the corresponding gene will be excluded from the gene-gene correlation estimation.
#' The default value is "all" (a special string which means no filtering).
#' @param nonnegative A logical variable. If TRUE, values < 0 in the synthetic data will be converted to 0. Default is TRUE (since the expression matrix is nonnegative).
#' @param nonzerovar A logical variable. If TRUE, for any gene with zero variance, a cell will be replaced with 1. This is designed for avoiding potential errors, for example, PCA. Default is FALSE.
#' @param sim_block_size Number of genes to simulate per bounded working block.
#' @param output_sparse If TRUE, return sparse simulated counts and avoid dense final materialization. The default follows whether the input assay is sparse.
#' @param sim_cell_block_size Number of cells per simulation block. For sparse atlas-scale inputs, parameter extraction and simulation are streamed by cell block to avoid full dense cell-by-gene parameter matrices.
#' @param return_model A logic variable. If TRUE, the marginal models and copula models will be returned. Default is FALSE.
#' @param simplify A logic variable. If TRUE, the fitted regression model will only keep the essential contains for \code{predict}, otherwise the fitted models can be VERY large. Default is FALSE.
#' @param parallelization A string indicating the specific parallelization function to use.
#' @param n_rep An integer number. The number of replicates of simulated new count matrix. Default is 1.
#' Must be one of 'mcmapply', 'bpmapply', or 'pbmcmapply', which corresponds to the parallelization function in the package
#' \code{parallel},\code{BiocParallel}, and \code{pbmcapply} respectively. The default value is 'mcmapply'.
#' @param BPPARAM A \code{MulticoreParam} object or NULL. When the parameter parallelization = 'mcmapply' or 'pbmcmapply',
#' this parameter must be NULL. When the parameter parallelization = 'bpmapply',  this parameter must be one of the
#' \code{MulticoreParam} object offered by the package 'BiocParallel. The default value is NULL.
#' @param trace A logic variable. If TRUE, the warning/error log and runtime for gam/gamlss
#' will be returned, FALSE otherwise. Default is FALSE.
#' @return A list with the components:
#' \describe{
#'   \item{\code{new_count}}{A matrix of the new simulated count (expression) matrix.}
#'   \item{\code{new_covariate}}{A data.frame of the covariates used for simulation. If \code{ncell} is unchanged, this contains the original covariates.}
#'   \item{\code{model_aic}}{The model AIC.}
#'   \item{\code{marginal_list}}{A list of marginal regression models if return_model = TRUE.}
#'   \item{\code{corr_list}}{A list of correlation models (conditional copulas) if return_model = TRUE.}
#' }
#' @examples
#' data(example_sce)
#' my_simu <- scdesign3(
#' sce = example_sce,
#' assay_use = "counts",
#' celltype = "cell_type",
#' pseudotime = "pseudotime",
#' spatial = NULL,
#' other_covariates = NULL,
#' mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
#' sigma_formula = "1",
#' family_use = "nb",
#' n_cores = 2,
#' usebam = FALSE,
#' edf_flexible = FALSE,
#' corr_formula = "pseudotime",
#' copula = "gaussian",
#' if_sparse = TRUE,
#' DT = TRUE,
#' pseudo_obs = FALSE,
#' ncell = 1000,
#' return_model = FALSE
#' )
#'
#' @export scdesign3
scdesign3 <- function(sce,
                      assay_use = "counts",
                      celltype,
                      pseudotime = NULL,
                      spatial = NULL,
                      other_covariates,
                      ncell = dim(sce)[2],
                      mu_formula,
                      sigma_formula = "1",
                      family_use = "nb",
                      n_cores = 2,
                      correlation_function = "default",
                      usebam = FALSE,
                      use_scglm = c("auto", "always", "never"),
                      scglm_method = c("auto", "categorical_closed_form", "categorical_irls", "irls", "newton_stein"),
                      scglm_batch_size = 256L,
                      edf_flexible = FALSE,
                      corr_formula,
                      empirical_quantile = FALSE,
                      copula = "gaussian",
                      if_sparse = FALSE,
                      fastmvn = FALSE,
                      DT = TRUE,
                      pseudo_obs = FALSE,
                      family_set = c("gauss", "indep"),
                      important_feature = "all",
                      nonnegative = TRUE,
                      nonzerovar = FALSE,
                      sim_block_size = 256L,
                      output_sparse = NULL,
                      sim_cell_block_size = getOption("scDesign3.sim_cell_block_size", 50000L),
                      return_model = FALSE,
                      simplify = FALSE,
                      parallelization = "mcmapply",
                      n_rep = 1,
                      BPPARAM = NULL,
                      trace = FALSE) {
  allowed_family <- c("binomial", "poisson", "nb", "zip", "zinb", "gaussian")
  allowed_copula <- c("gaussian", "vine")
  allowed_parallelization <- c("mcmapply", "bpmapply", "pbmcmapply")
  allowed_correlation_function <- c("default", "coop")
  use_scglm <- match.arg(use_scglm)
  scglm_method <- match.arg(scglm_method)

  if (!copula %in% allowed_copula) {
    stop("copula must be one of 'gaussian' or 'vine'.")
  }
  if (!parallelization %in% allowed_parallelization) {
    stop("parallelization must be one of 'mcmapply', 'bpmapply', or 'pbmcmapply'.")
  }
  if (!correlation_function %in% allowed_correlation_function) {
    stop("correlation_function must be one of 'default' or 'coop'.")
  }
  if (length(family_use) < 1 || any(!family_use %in% allowed_family)) {
    stop("family_use must contain only 'binomial', 'poisson', 'nb', 'zip', 'zinb', or 'gaussian'.")
  }
  if (n_rep < 1 || length(n_rep) != 1 || !is.numeric(n_rep)) {
    stop("n_rep must be a positive integer.")
  }
  if (n_rep != as.integer(n_rep)) {
    stop("n_rep must be an integer.")
  }
  if (empirical_quantile && ncell != dim(sce)[2]) {
    stop("empirical_quantile = TRUE only works when ncell equals the number of cells in sce.")
  }

  message("Input Data Construction Start")

  input_data <- construct_data(
    sce = sce,
    assay_use = assay_use,
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    ncell = ncell,
    corr_by = corr_formula,
    parallelization = parallelization,
    BPPARAM = BPPARAM
  )
  message("Input Data Construction End")

  message("Start Marginal Fitting")
  marginal_res <- fit_marginal(
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    n_cores = n_cores,
    data = input_data,
    family_use = family_use,
    usebam = usebam,
    use_scglm = use_scglm,
    scglm_method = scglm_method,
    scglm_batch_size = scglm_batch_size,
    edf_flexible = edf_flexible,
    parallelization = parallelization,
    BPPARAM = BPPARAM,
    trace = trace, 
    simplify = simplify
  )
  message("Marginal Fitting End")

  if(empirical_quantile == TRUE) {
    message("Extract Empirical Quantile Matrices")
    copula_res <- fit_copula(
      sce = sce,
      assay_use = assay_use,
      input_data = input_data$dat,
      marginal_list = marginal_res,
      family_use = family_use,
      empirical_quantile = TRUE,
      copula = copula,
      DT = DT,
      pseudo_obs = pseudo_obs,
      family_set = family_set,
      n_cores = n_cores,
      important_feature = important_feature,
      if_sparse = if_sparse,
      parallelization = parallelization,
      BPPARAM = BPPARAM
    )
  } else {
    message("Start Copula Fitting")
    copula_res <- fit_copula(
      sce = sce,
      assay_use = assay_use,
      input_data = input_data$dat,
      marginal_list = marginal_res,
      family_use = family_use,
      copula = copula,
      DT = DT,
      pseudo_obs = pseudo_obs,
      family_set = family_set,
      n_cores = n_cores,
      correlation_function = correlation_function,
      important_feature = important_feature,
      if_sparse = if_sparse,
      parallelization = parallelization,
      BPPARAM = BPPARAM
    )
    message("Copula Fitting End")
  }
  
  

  stream_simulation <- .scdesign3_should_stream_simulation(
    sce = sce,
    assay_use = assay_use,
    new_covariate = input_data$new_covariate,
    output_sparse = output_sparse,
    sim_cell_block_size = sim_cell_block_size
  )

  if (stream_simulation) {
    message("Start Blocked Parameter Extraction and Data Generation")
    new_count <- .scdesign3_blocked_extract_and_simulate(
      sce = sce,
      assay_use = assay_use,
      marginal_res = marginal_res,
      copula_res = copula_res,
      input_data = input_data,
      family_use = family_use,
      n_cores = n_cores,
      fastmvn = fastmvn,
      nonnegative = nonnegative,
      nonzerovar = nonzerovar,
      empirical_quantile = empirical_quantile,
      parallelization = parallelization,
      BPPARAM = BPPARAM,
      sim_block_size = sim_block_size,
      output_sparse = output_sparse,
      sim_cell_block_size = sim_cell_block_size,
      n_rep = n_rep
    )
    message("Blocked Parameter Extraction and Data Generation End")
  } else {
    message("Start Parameter Extraction")
    para_list <- extract_para(
      sce = sce,
      assay_use = assay_use,
      marginal_list = marginal_res,
      n_cores = n_cores,
      family_use = family_use,
      new_covariate = input_data$new_covariate,
      parallelization = parallelization,
      BPPARAM = BPPARAM,
      data = input_data$dat
    )
    message("Parameter Extraction End")

    message("Start Generate New Data")
    new_count <- .scdesign3_simulate_from_para(
      sce = sce,
      assay_use = assay_use,
      para_list = para_list,
      copula_res = copula_res,
      input_data = input_data,
      family_use = family_use,
      n_cores = n_cores,
      fastmvn = fastmvn,
      nonnegative = nonnegative,
      nonzerovar = nonzerovar,
      empirical_quantile = empirical_quantile,
      parallelization = parallelization,
      BPPARAM = BPPARAM,
      sim_block_size = sim_block_size,
      output_sparse = output_sparse,
      n_rep = n_rep
    )
  }
  
  message("New Data Generating End")

  scdesign3_res <- list(
    new_count = new_count,
    new_covariate = input_data$new_covariate,
    model_aic = copula_res$model_aic,
    model_bic = copula_res$model_bic,
    marginal_list = if (return_model)
      marginal_res
    else
      NULL,
    corr_list = if (return_model)
      copula_res$copula_list
    else
      NULL
  )
  return(scdesign3_res)
}

.scdesign3_simulate_from_para <- function(sce,
                                          assay_use,
                                          para_list,
                                          copula_res,
                                          input_data,
                                          family_use,
                                          n_cores,
                                          fastmvn,
                                          nonnegative,
                                          nonzerovar,
                                          empirical_quantile,
                                          parallelization,
                                          BPPARAM,
                                          sim_block_size,
                                          output_sparse,
                                          n_rep) {
  simulate_one <- function(quantile_mat, copula_list) {
    simu_new(
      sce = sce,
      assay_use= assay_use,
      mean_mat = para_list$mean_mat,
      sigma_mat = para_list$sigma_mat,
      zero_mat = para_list$zero_mat,
      quantile_mat = quantile_mat,
      copula_list = copula_list,
      n_cores = n_cores,
      fastmvn = fastmvn,
      family_use = family_use,
      nonnegative = nonnegative,
      nonzerovar = nonzerovar,
      input_data = input_data$dat,
      new_covariate = input_data$new_covariate,
      important_feature = copula_res$important_feature,
      parallelization = parallelization,
      BPPARAM = BPPARAM,
      filtered_gene = input_data$filtered_gene,
      sim_block_size = sim_block_size,
      output_sparse = output_sparse
    )
  }

  if (empirical_quantile) {
    return(simulate_one(copula_res$quantile_mat, NULL))
  }
  if (n_rep == 1) {
    return(simulate_one(NULL, copula_res$copula_list))
  }
  lapply(seq_len(n_rep), function(i) simulate_one(NULL, copula_res$copula_list))
}

.scdesign3_should_stream_simulation <- function(sce,
                                                assay_use,
                                                new_covariate,
                                                output_sparse,
                                                sim_cell_block_size) {
  sim_cell_block_size <- suppressWarnings(as.integer(sim_cell_block_size[1L]))
  if (!is.finite(sim_cell_block_size) || sim_cell_block_size < 1L) {
    return(FALSE)
  }
  if (is.null(output_sparse)) {
    output_sparse <- methods::is(SummarizedExperiment::assay(sce, assay_use), "sparseMatrix")
  }
  if (!isTRUE(output_sparse)) {
    return(FALSE)
  }
  total_cells <- if (is.null(new_covariate)) {
    dim(sce)[2]
  } else {
    nrow(new_covariate)
  }
  total_cells > sim_cell_block_size
}

.scdesign3_blocked_extract_and_simulate <- function(sce,
                                                   assay_use,
                                                   marginal_res,
                                                   copula_res,
                                                   input_data,
                                                   family_use,
                                                   n_cores,
                                                   fastmvn,
                                                   nonnegative,
                                                   nonzerovar,
                                                   empirical_quantile,
                                                   parallelization,
                                                   BPPARAM,
                                                   sim_block_size,
                                                   output_sparse,
                                                   sim_cell_block_size,
                                                   n_rep) {
  new_covariate_all <- input_data$new_covariate
  if (is.null(new_covariate_all)) {
    new_covariate_all <- input_data$dat
  }
  if (!"corr_group" %in% colnames(new_covariate_all) &&
      "corr_group" %in% colnames(input_data$dat)) {
    corr_idx <- match(rownames(new_covariate_all), rownames(input_data$dat))
    new_covariate_all$corr_group <- input_data$dat$corr_group[corr_idx]
  }
  total_cells <- nrow(new_covariate_all)
  blocks <- split(seq_len(total_cells), ceiling(seq_len(total_cells) / sim_cell_block_size))

  simulate_rep <- function(rep_id) {
    block_counts <- lapply(blocks, function(idx) {
      new_covariate_block <- new_covariate_all[idx, , drop = FALSE]
      para_block <- extract_para(
        sce = sce,
        assay_use = assay_use,
        marginal_list = marginal_res,
        n_cores = n_cores,
        family_use = family_use,
        new_covariate = new_covariate_block,
        parallelization = parallelization,
        BPPARAM = BPPARAM,
        data = input_data$dat
      )
      quantile_block <- if (empirical_quantile) {
        copula_res$quantile_mat[idx, , drop = FALSE]
      } else {
        NULL
      }
      simu_new(
        sce = sce,
        assay_use= assay_use,
        mean_mat = para_block$mean_mat,
        sigma_mat = para_block$sigma_mat,
        zero_mat = para_block$zero_mat,
        quantile_mat = quantile_block,
        copula_list = if (empirical_quantile) NULL else copula_res$copula_list,
        n_cores = n_cores,
        fastmvn = fastmvn,
        family_use = family_use,
        nonnegative = nonnegative,
        nonzerovar = nonzerovar,
        input_data = input_data$dat,
        new_covariate = new_covariate_block,
        important_feature = copula_res$important_feature,
        parallelization = parallelization,
        BPPARAM = BPPARAM,
        filtered_gene = input_data$filtered_gene,
        sim_block_size = sim_block_size,
        output_sparse = output_sparse
      )
    })
    do.call(cbind, block_counts)
  }

  if (n_rep == 1) {
    return(simulate_rep(1L))
  }
  lapply(seq_len(n_rep), simulate_rep)
}
