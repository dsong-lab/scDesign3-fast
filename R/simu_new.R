#' Simulate new data
#'
#' \code{simu_new} generates new simulated data based on fitted marginal and copula models.
#'
#' The function takes the new covariate (if use) from \code{\link{construct_data}},
#' parameter matrices from \code{\link{extract_para}} and multivariate Unifs from \code{\link{fit_copula}}.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assay_use A string which indicates the assay you will use in the sce. Default is 'counts'.
#' @param mean_mat A cell by feature matrix of the mean parameter.
#' @param sigma_mat A cell by feature matrix of the sigma parameter.
#' @param zero_mat A cell by feature matrix of the zero-inflation parameter.
#' @param quantile_mat A cell by feature matrix of the multivariate quantile.
#' @param copula_list A list of copulas for generating the multivariate quantile matrix. If provided, the \code{quantile_mat} must be NULL.
#' @param n_cores An integer. The number of cores to use.
#' @param fastmvn An logical variable. If TRUE, the sampling of multivariate Gaussian is done by \code{mvnfast}, otherwise by \code{mvtnorm}. Default is FALSE.
#' @param family_use A string of the marginal distribution.
#' Must be one of 'poisson', "binomial", 'nb', 'zip', 'zinb' or 'gaussian'.
#' @param nonnegative A logical variable. If TRUE, values < 0 in the synthetic data will be converted to 0. Default is TRUE (since the expression matrix is nonnegative).
#' @param nonzerovar A logical variable. If TRUE, for any gene with zero variance, a cell will be replaced with 1. This is designed for avoiding potential errors, for example, PCA.
#' @param input_data A input count matrix.
#' @param new_covariate A data.frame which contains covariates of targeted simulated data from  \code{\link{construct_data}}.
#' @param important_feature important_feature A string or vector which indicates whether a gene will be used in correlation estimation or not. If this is a string, then
#' this string must be either "all" (using all genes) or "auto", which indicates that the genes will be automatically selected based on the proportion of zero expression across cells
#' for each gene. Gene with zero proportion greater than 0.8 will be excluded form gene-gene correlation estimation. If this is a vector, then this should
#' be a logical vector with length equal to the number of genes in \code{sce}. \code{TRUE} in the logical vector means the corresponding gene will be included in
#' gene-gene correlation estimation and \code{FALSE} in the logical vector means the corresponding gene will be excluded from the gene-gene correlation estimation.
#' The default value for is "all".
#' @param parallelization A string indicating the specific parallelization function to use.
#' Must be one of 'mcmapply', 'bpmapply', or 'pbmcmapply', which corresponds to the parallelization function in the package
#' \code{parallel},\code{BiocParallel}, and \code{pbmcapply} respectively. The default value is 'mcmapply'.
#' @param BPPARAM A \code{MulticoreParam} object or NULL. When the parameter parallelization = 'mcmapply' or 'pbmcmapply',
#' this parameter must be NULL. When the parameter parallelization = 'bpmapply',  this parameter must be one of the
#' \code{MulticoreParam} object offered by the package 'BiocParallel. The default value is NULL.
#' @param filtered_gene A vector or NULL which contains genes that are excluded in the marginal and copula fitting 
#' steps because these genes only express in less than two cells. This can be obtain from  \code{\link{construct_data}}
#' @param sim_block_size Number of genes to simulate per bounded working block.
#' @param output_sparse If TRUE, return a sparse feature-by-cell matrix and avoid dense final materialization. The default follows whether the input assay is sparse.
#' @return A feature by cell matrix of the new simulated count (expression) matrix or sparse matrix.
#' @examples
#'   data(example_sce)
#'   my_data <- construct_data(
#'   sce = example_sce,
#'   assay_use = "counts",
#'   celltype = "cell_type",
#'   pseudotime = "pseudotime",
#'   spatial = NULL,
#'   other_covariates = NULL,
#'   corr_by = "1"
#'   )
#'   my_marginal <- fit_marginal(
#'   data = my_data,
#'   mu_formula = "s(pseudotime, bs = 'cr', k = 10)",
#'   sigma_formula = "1",
#'   family_use = "nb",
#'   n_cores = 1,
#'   usebam = FALSE
#'   )
#'   my_copula <- fit_copula(
#'   sce = example_sce,
#'   assay_use = "counts",
#'   marginal_list = my_marginal,
#'   family_use = c(rep("nb", 5), rep("zip", 5)),
#'   copula = "vine",
#'   n_cores = 1,
#'   input_data = my_data$dat
#'   )
#'   my_para <- extract_para(
#'     sce = example_sce,
#'     marginal_list = my_marginal,
#'     n_cores = 1,
#'     family_use = c(rep("nb", 5), rep("zip", 5)),
#'     new_covariate = my_data$new_covariate,
#'     data = my_data$dat
#'   )
#'   my_newcount <- simu_new(
#'   sce = example_sce,
#'   mean_mat = my_para$mean_mat,
#'   sigma_mat = my_para$sigma_mat,
#'   zero_mat = my_para$zero_mat,
#'   quantile_mat = NULL,
#'   copula_list = my_copula$copula_list,
#'   n_cores = 1,
#'   family_use = c(rep("nb", 5), rep("zip", 5)),
#'   input_data = my_data$dat,
#'   new_covariate = my_data$new_covariate,
#'   important_feature = my_copula$important_feature,
#'   filtered_gene = my_data$filtered_gene
#'   )
#'
#' @export simu_new

simu_new <- function(sce,
                     assay_use = "counts",
                     mean_mat,
                     sigma_mat,
                     zero_mat,
                     quantile_mat = NULL,
                     copula_list,
                     n_cores,
                     fastmvn = FALSE,
                     family_use,
                     nonnegative = TRUE,
                     nonzerovar = FALSE,
                     input_data,
                     new_covariate,
                     important_feature = "all",
                     parallelization = "mcmapply",
                     BPPARAM = NULL,
                     filtered_gene,
                     sim_block_size = 256L,
                     output_sparse = NULL){
  if(!is.null(quantile_mat) & !is.null(copula_list)) {
    stop("You can only provide either the quantile_mat or the copula_list!")
  }

  # check if user inputted new covariates
  data_temp <- input_data[,colnames(new_covariate), drop = FALSE]
  if(identical(data_temp, new_covariate)){
    new_covariate <- NULL
  }
  
  qc_gene_idx <- which(!rownames(sce) %in% filtered_gene)
  qc_gene_names <- rownames(sce)[qc_gene_idx]
  if(length(family_use) != 1){
    family_use <- family_use[qc_gene_idx]
  } else {
    family_use <- rep(family_use, length(qc_gene_idx))
  }

  if(is.null(new_covariate)){
    total_cells <- dim(sce)[2]
    cell_names <- colnames(sce)
  }else{
    total_cells <- dim(new_covariate)[1]
    cell_names <- rownames(new_covariate)
  }
  if (is.null(output_sparse)) {
    output_sparse <- methods::is(SummarizedExperiment::assay(sce, assay_use), "sparseMatrix")
  } else {
    output_sparse <- isTRUE(output_sparse)
  }

  if(!is.null(quantile_mat)) {
    message("Multivariate quantile matrix is provided")
    quantile_mat <- .scdesign3_align_quantile_matrix(
      quantile_mat = quantile_mat,
      cell_names = cell_names,
      gene_names = qc_gene_names
    )
  } else {
    message("Use Copula to sample a multivariate quantile matrix")
      
    group_index <- unique(input_data$corr_group)
    corr_group <- as.data.frame(input_data$corr_group)
    colnames(corr_group) <- "corr_group"
    ngene <- length(qc_gene_idx)
    if (is.null(new_covariate)) {
      new_corr_group <- NULL
    } else{
      new_corr_group <- as.data.frame(new_covariate$corr_group)
      colnames(new_corr_group) <- "corr_group"
    }
    ind <- group_index[1] == "ind"
    newmvn.list <-
      lapply(group_index, function(x,
                                   sce,
                                   corr_group,
                                   new_corr_group,
                                   ind,
                                   n_cores,
                                   copula_list) {
        message(paste0("Sample Copula group ", x, " starts"))
        curr_index <- which(corr_group[, 1] == x)
        if (is.null(new_covariate)) {
          curr_ncell <- length(curr_index)
          curr_ncell_idx <- curr_index
        } else{
          curr_ncell <- length(which(new_corr_group[, 1] == x))
          curr_ncell_idx <-which(new_corr_group[, 1] == x)
          #paste0("Cell", which(new_corr_group[, 1] == x))
        }
        cor.mat <- copula_list[[x]]
        
        if(curr_ncell == 0) {
          new_mvu <- NULL
        } else {
          if (inherits(cor.mat, "scdesign3_gaussian_fast_copula")) {
            new_mvu <- sampleMVN(
              n = curr_ncell,
              Sigma = cor.mat,
              n_cores = n_cores,
              fastmvn = fastmvn
            )
            new_mvu <- as.matrix(new_mvu)
            new_mvu <- new_mvu[, qc_gene_names, drop = FALSE]
            rownames(new_mvu) <- curr_ncell_idx
          } else if (methods::is(cor.mat, "matrix") | methods::is(cor.mat, "dsCMatrix")) {
            #message(paste0("Group ", group_index, " Start"))
            
            #message("Sample MVN")
            #sample from mvn for important genes only
            cache_names <- attr(cor.mat, "scdesign3_gaussian_factor_names", exact = TRUE)
            if (!is.null(cache_names)) {
              corr_gene <- cache_names
            } else {
              corr_gene <- .scdesign3_correlated_gene_names(cor.mat)
            }
            corr_gene_idx <- colnames(cor.mat) %in% corr_gene
            if(length(corr_gene)!=0) {
              sigma_sample <- if (!is.null(cache_names)) {
                cor.mat
              } else {
                cor.mat[corr_gene, corr_gene, drop = FALSE]
              }
              new_mvn_important <- sampleMVN(
                n = curr_ncell,
                Sigma = sigma_sample,
                n_cores = n_cores,
                fastmvn = fastmvn
              )
              if (is.null(dim(new_mvn_important))) {
                new_mvn_important <- matrix(new_mvn_important, ncol = length(corr_gene))
              } else {
                new_mvn_important <- as.matrix(new_mvn_important)
              }
            colnames(new_mvn_important) <- corr_gene} else {
              new_mvn_important <- NULL
            }
            #message("MVN Sampling End")
            ind_gene <- colnames(cor.mat)[which(corr_gene_idx==FALSE)]
            if(length(ind_gene) > 0){
              new_mvn_non_important_mat <- matrix(
                stats::rnorm(curr_ncell * length(ind_gene)),
                nrow = curr_ncell,
                ncol = length(ind_gene)
              )
              colnames(new_mvn_non_important_mat) <- ind_gene
              
              mvnrvq <- matrix(
                stats::pnorm(new_mvn_non_important_mat),
                nrow = nrow(new_mvn_non_important_mat),
                ncol = ncol(new_mvn_non_important_mat),
                dimnames = dimnames(new_mvn_non_important_mat)
              )
              new_mvu <- cbind(new_mvn_important, mvnrvq)
              new_mvu <- new_mvu[,colnames(cor.mat)]
            }else{
              new_mvu <- new_mvn_important
            }

            rownames(new_mvu) <- curr_ncell_idx
          } else if (methods::is(cor.mat, "vinecop")) {
            new_mvu <- matrix(0, nrow = curr_ncell, ncol = ngene)
            #message("Sampling Vine Copula Starts")
            mvu <- rvinecopulib::rvinecop(
              curr_ncell,
              vine = cor.mat,
              cores = n_cores,
              qrng = TRUE
            )
            new_mvu[, which(important_feature)] <- mvu
            if(length(which(important_feature)) != ngene){
              cor.mat <- diag(rep(1, length(which(!important_feature))))
              mvu2 <- sampleMVN(n = curr_ncell,
                                Sigma = cor.mat,
                                n_cores = n_cores,
                                fastmvn = fastmvn)
              new_mvu[, which(!important_feature)] <- mvu2
            }
            #message("Sampling Vine Copula Ends")
            rownames(new_mvu) <- curr_ncell_idx
          } else if (ind) {
            "Use independent copula (random Unif)."
            new_mvu <-
              matrix(data = stats::runif(curr_ncell * ngene),
                     nrow = curr_ncell)
            rownames(new_mvu) <- curr_ncell_idx
          } else{
            stop("Copula must be one from 'vine' or 'gaussian', or assume gene-gene is independent")
          }
        }
        return(
          list(
            new_mvu = new_mvu
          )
        )
      }, sce = sce, ind = ind, n_cores = n_cores, corr_group = corr_group, new_corr_group = new_corr_group, copula_list = copula_list)
    
    newmvn <-
      do.call(rbind, lapply(newmvn.list, function(x)
        x$new_mvu))
    newmvn[as.numeric(rownames(newmvn)),] <- newmvn
    rownames(newmvn) <- as.character(seq_len(dim(newmvn)[1]))
    colnames(newmvn) <- qc_gene_names
    quantile_mat <- newmvn
  }

  ## New count
  new_count <- .scdesign3_quantile_to_counts(
    mean_mat = mean_mat,
    sigma_mat = sigma_mat,
    zero_mat = zero_mat,
    quantile_mat = quantile_mat,
    gene_idx = qc_gene_idx,
    family_use = family_use,
    cell_names = cell_names,
    n_cores = n_cores,
    parallelization = parallelization,
    BPPARAM = BPPARAM,
    block_size = sim_block_size,
    output_sparse = output_sparse
  )

  if(length(qc_gene_idx) < dim(sce)[1]){
    if (isTRUE(output_sparse)) {
      temp_count <- Matrix::Matrix(0, nrow = dim(sce)[1], ncol = total_cells, sparse = TRUE)
    } else {
      temp_count <- matrix(0, nrow = dim(sce)[1], ncol = total_cells)
    }
    rownames(temp_count) <- rownames(sce)
    colnames(temp_count) <- cell_names
    temp_count[rownames(new_count), colnames(new_count)] <- new_count
    new_count <- temp_count
  }

  if(nonnegative) {
    if (methods::is(new_count, "sparseMatrix")) {
      new_count@x[new_count@x < 0] <- 0
      new_count <- Matrix::drop0(new_count)
    } else {
      new_count[new_count < 0] <- 0
    }
  }

  if(nonzerovar) {
    new_count <- .scdesign3_fix_zero_variance_rows(new_count, qc_gene_idx)
  }

  return(new_count)
}

.scdesign3_align_quantile_matrix <- function(quantile_mat, cell_names, gene_names) {
  quantile_mat <- as.matrix(quantile_mat)
  out <- matrix(
    NA_real_,
    nrow = length(cell_names),
    ncol = length(gene_names),
    dimnames = list(cell_names, gene_names)
  )

  if (is.null(rownames(quantile_mat))) {
    row_to <- seq_len(min(nrow(out), nrow(quantile_mat)))
    row_from <- row_to
  } else {
    row_match <- match(rownames(quantile_mat), cell_names)
    if (all(is.na(row_match))) {
      row_to <- seq_len(min(nrow(out), nrow(quantile_mat)))
      row_from <- row_to
    } else {
      keep <- !is.na(row_match)
      row_to <- row_match[keep]
      row_from <- which(keep)
    }
  }

  if (is.null(colnames(quantile_mat))) {
    col_to <- seq_len(min(ncol(out), ncol(quantile_mat)))
    col_from <- col_to
  } else {
    col_match <- match(colnames(quantile_mat), gene_names)
    if (all(is.na(col_match))) {
      col_to <- seq_len(min(ncol(out), ncol(quantile_mat)))
      col_from <- col_to
    } else {
      keep <- !is.na(col_match)
      col_to <- col_match[keep]
      col_from <- which(keep)
    }
  }

  if (length(row_to) > 0L && length(col_to) > 0L) {
    out[row_to, col_to] <- quantile_mat[row_from, col_from, drop = FALSE]
  }
  out
}

.scdesign3_quantile_to_counts <- function(mean_mat,
                                          sigma_mat,
                                          zero_mat,
                                          quantile_mat,
                                          gene_idx,
                                          family_use,
                                          cell_names,
                                          n_cores,
                                          parallelization,
                                          BPPARAM,
                                          block_size = 256L,
                                          output_sparse = FALSE) {
  block_size <- as.integer(block_size[1L])
  if (!is.finite(block_size) || block_size < 1L) {
    block_size <- 256L
  }
  n_cores <- as.integer(n_cores[1L])
  if (!is.finite(n_cores) || n_cores < 1L) {
    n_cores <- 1L
  }
  if (is.null(parallelization) || !nzchar(parallelization[1L])) {
    parallelization <- "mcmapply"
  }
  if (length(gene_idx) == 0L) {
    if (isTRUE(output_sparse)) {
      return(Matrix::Matrix(0, nrow = 0L, ncol = length(cell_names), sparse = TRUE))
    }
    return(matrix(0, nrow = 0L, ncol = length(cell_names), dimnames = list(NULL, cell_names)))
  }
  gene_pos <- seq_along(gene_idx)
  blocks <- split(gene_pos, ceiling(gene_pos / block_size))

  block_fun <- function(pos) {
    genes <- gene_idx[pos]
    q <- as.matrix(quantile_mat[, pos, drop = FALSE])
    mu <- as.matrix(mean_mat[, genes, drop = FALSE])
    fam <- family_use[pos]

    out <- matrix(0, nrow = nrow(q), ncol = ncol(q))
    for (family_i in unique(fam)) {
      cols <- which(fam == family_i)
      sigma <- if (family_i %in% c("gaussian", "nb", "zinb")) {
        as.matrix(sigma_mat[, genes[cols], drop = FALSE])
      } else {
        NULL
      }
      zero <- if (family_i %in% c("zip", "zinb")) {
        as.matrix(zero_mat[, genes[cols], drop = FALSE])
      } else {
        NULL
      }
      out[, cols] <- .scdesign3_qfamily_matrix(
        y = family_i,
        p = q[, cols, drop = FALSE],
        mu = mu[, cols, drop = FALSE],
        sigma = sigma,
        zero = zero
      )
    }
    out
  }

  if (.Platform$OS.type == "windows") {
    BPPARAM <- BiocParallel::SnowParam()
    parallelization <- "bpmapply"
  }

  if (parallelization != "bpmapply" && !(parallelization == "pbmcmapply" && n_cores > 1L) && n_cores <= 1L) {
    if (isTRUE(output_sparse)) {
      out_blocks <- lapply(blocks, function(pos) {
        block <- Matrix::Matrix(t(block_fun(pos)), sparse = TRUE)
        rownames(block) <- colnames(mean_mat)[gene_idx[pos]]
        colnames(block) <- cell_names
        block
      })
      return(do.call(rbind, out_blocks))
    } else {
      out <- matrix(
        0,
        nrow = length(gene_idx),
        ncol = nrow(quantile_mat),
        dimnames = list(colnames(mean_mat)[gene_idx], cell_names)
      )
      for (pos in blocks) {
        out[pos, ] <- t(block_fun(pos))
      }
      return(out)
    }
  }

  block_results <- if (parallelization == "bpmapply") {
    if (is.null(BPPARAM)) {
      BPPARAM <- BiocParallel::SerialParam()
    }
    if (!methods::is(BPPARAM, "SerialParam")) {
      BPPARAM$workers <- n_cores
    }
    BiocParallel::bplapply(blocks, block_fun, BPPARAM = BPPARAM)
  } else if (parallelization == "pbmcmapply" && requireNamespace("pbmcapply", quietly = TRUE)) {
    pbmcapply::pbmclapply(blocks, block_fun, mc.cores = n_cores)
  } else if (n_cores > 1L) {
    parallel::mclapply(blocks, block_fun, mc.cores = n_cores)
  } else {
    lapply(blocks, block_fun)
  }

  if (isTRUE(output_sparse)) {
    sparse_blocks <- lapply(seq_along(blocks), function(i) {
      pos <- blocks[[i]]
      block <- Matrix::Matrix(t(block_results[[i]]), sparse = TRUE)
      rownames(block) <- colnames(mean_mat)[gene_idx[pos]]
      colnames(block) <- cell_names
      block
    })
    return(do.call(rbind, sparse_blocks))
  }
  out <- do.call(rbind, lapply(block_results, t))
  rownames(out) <- colnames(mean_mat)[gene_idx]
  colnames(out) <- cell_names
  out
}

.scdesign3_fix_zero_variance_rows <- function(new_count, qc_gene_idx) {
  if (methods::is(new_count, "sparseMatrix")) {
    row_means <- Matrix::rowMeans(new_count[qc_gene_idx, , drop = FALSE])
    row_sq_means <- Matrix::rowMeans(new_count[qc_gene_idx, , drop = FALSE]^2)
    row_vars <- as.numeric(row_sq_means - row_means^2)
  } else {
    row_vars <- matrixStats::rowVars(new_count[qc_gene_idx, , drop = FALSE])
  }
  if(sum(row_vars == 0) > 0) {
    message("Some genes have zero variance. Replace a random one with 1.")
    row_vars_index <- qc_gene_idx[which(row_vars == 0)]
    col_index <- seq_len(dim(new_count)[2])
    for(i in row_vars_index) {
      new_count[i, sample(col_index, 1)] <- 1
    }
  }
  if (methods::is(new_count, "sparseMatrix")) {
    new_count <- Matrix::drop0(new_count)
  }
  new_count
}

.scdesign3_qfamily_matrix <- function(y, p, mu, sigma, zero) {
  nr <- nrow(p)
  nc <- ncol(p)
  p_vec <- as.vector(p)
  mu_vec <- as.vector(mu)

  if (y == "binomial") {
    out <- stats::qbinom(p = p_vec, prob = mu_vec, size = 1)
  } else if (y == "poisson") {
    out <- stats::qpois(p = p_vec, lambda = mu_vec)
  } else if (y == "gaussian") {
    sigma_vec <- as.vector(sigma)
    out <- gamlss.dist::qNO(p = p_vec, mu = mu_vec, sigma = abs(sigma_vec))
  } else if (y == "nb") {
    sigma_vec <- as.vector(sigma)
    zero_mu <- is.na(mu_vec) | mu_vec <= 0
    mu_safe <- mu_vec
    sigma_safe <- sigma_vec
    mu_safe[zero_mu] <- .Machine$double.eps
    sigma_safe[!is.finite(sigma_safe) | sigma_safe <= 0] <- .Machine$double.eps
    out <- gamlss.dist::qNBI(p = p_vec, mu = mu_safe, sigma = sigma_safe)
    out[zero_mu] <- 0
  } else if (y == "zip") {
    zero_vec <- as.vector(zero)
    out <- gamlss.dist::qZIP(
      p = p_vec,
      mu = mu_vec,
      sigma = ifelse(zero_vec != 0, zero_vec, 2.2e-16)
    )
  } else if (y == "zinb") {
    sigma_vec <- as.vector(sigma)
    zero_vec <- as.vector(zero)
    out <- gamlss.dist::qZINBI(
      p = p_vec,
      mu = mu_vec,
      sigma = sigma_vec,
      nu = ifelse(zero_vec != 0, zero_vec, 2.2e-16)
    )
  } else {
    stop("Distribution of gamlss must be one of gaussian, poisson, nb, zip or zinb!")
  }

  out <- matrix(out, nrow = nr, ncol = nc)
  out[is.na(mu) | mu == 0] <- 0
  out
}
