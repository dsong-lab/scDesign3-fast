test_that("vectorized marginal CDF transform matches legacy row-apply logic", {
  legacy_uniform <- function(y, observed, mean, theta, zero = 0, DT = TRUE) {
    family_frame <- cbind(observed, mean, theta, zero)
    if (y == "poisson") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::ppois(x[1], lambda = x[2])
      })
    } else if (y == "nb") {
      pvec <- apply(family_frame, 1, function(x) {
        stats::pnbinom(x[1], mu = x[2], size = x[3])
      })
    } else {
      stop("unsupported test family")
    }

    if (DT) {
      if (y == "poisson") {
        pvec2 <- apply(family_frame, 1, function(x) {
          stats::ppois(x[1] - 1, lambda = x[2]) * as.integer(x[1] > 0)
        })
      } else {
        pvec2 <- apply(family_frame, 1, function(x) {
          stats::pnbinom(x[1] - 1, mu = x[2], size = x[3]) * as.integer(x[1] > 0)
        })
      }
      v <- stats::runif(length(mean))
      pvec <- pvec * v + pvec2 * (1 - v)
    }
    names(pvec) <- names(mean)
    pvec
  }

  observed <- c(0, 1, 4, 2, 7)
  mean <- stats::setNames(c(0.4, 1.2, 3.5, 2.1, 6.8), paste0("cell", seq_along(observed)))
  theta <- c(2.5, 2.5, 3, 4, 5)

  set.seed(401)
  old_poisson <- legacy_uniform("poisson", observed, mean, theta, DT = TRUE)
  set.seed(401)
  new_poisson <- .scdesign3_marginal_uniform("poisson", observed, mean, theta, DT = TRUE)

  set.seed(402)
  old_nb <- legacy_uniform("nb", observed, mean, theta, DT = TRUE)
  set.seed(402)
  new_nb <- .scdesign3_marginal_uniform("nb", observed, mean, theta, DT = TRUE)

  expect_equal(new_poisson, old_poisson, tolerance = 1e-15)
  expect_equal(new_nb, old_nb, tolerance = 1e-15)
})

test_that("gene zero proportion helper is sparse-safe", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Matrix")

  mat_dense <- matrix(
    c(0, 1, 0, 0, 2, 3, 0, 0, 0, 4, 0, 5),
    nrow = 3,
    byrow = TRUE
  )
  rownames(mat_dense) <- paste0("gene", seq_len(nrow(mat_dense)))
  colnames(mat_dense) <- paste0("cell", seq_len(ncol(mat_dense)))
  sce_dense <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat_dense))
  sce_sparse <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(mat_dense, sparse = TRUE))
  )

  expect_equal(
    .scdesign3_gene_zero_prop(sce_sparse, "counts"),
    .scdesign3_gene_zero_prop(sce_dense, "counts")
  )
})
