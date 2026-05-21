#' @importFrom utils getFromNamespace
NULL

.scdesign3_fit_marginal_scglm <- function(count_mat,
                                          dat_cov,
                                          filtered_gene,
                                          feature_names,
                                          family_use,
                                          mu_formula,
                                          sigma_formula,
                                          predictor,
                                          use_scglm,
                                          scglm_method,
                                          scglm_batch_size,
                                          n_cores,
                                          trace,
                                          filter_cells) {
  if (identical(use_scglm, "never")) {
    return(NULL)
  }

  unsupported <- function(reason) {
    if (identical(use_scglm, "always")) {
      stop(reason, call. = FALSE)
    }
    NULL
  }

  if (isTRUE(filter_cells)) {
    return(unsupported("The scGLM backend does not support gene-specific `filter_cells`; use the default backend."))
  }
  if (length(unique(family_use)) != 1L) {
    return(unsupported("The scGLM backend currently requires one shared `family_use` for all features."))
  }

  family_key <- unique(family_use)
  categorical_fit <- .scdesign3_fit_categorical_gamlss(
    count_mat = count_mat,
    dat_cov = dat_cov,
    filtered_gene = filtered_gene,
    feature_names = feature_names,
    family_key = family_key,
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    predictor = predictor,
    use_scglm = use_scglm,
    trace = trace
  )
  if (!is.null(categorical_fit)) {
    return(categorical_fit)
  }

  smooth_fit <- .scdesign3_fit_shared_fixed_smooth(
    count_mat = count_mat,
    dat_cov = dat_cov,
    filtered_gene = filtered_gene,
    feature_names = feature_names,
    family_key = family_key,
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    predictor = predictor,
    use_scglm = use_scglm,
    scglm_method = scglm_method,
    scglm_batch_size = scglm_batch_size,
    n_cores = n_cores,
    trace = trace
  )
  if (!is.null(smooth_fit)) {
    return(smooth_fit)
  }
  if (.scdesign3_formula_has_smooth(mu_formula) ||
      .scdesign3_formula_has_smooth(sigma_formula)) {
    return(unsupported(
      "The scGLM smooth backend is used only for fixed-basis or fixed-penalty smooths with `sigma_formula = \"1\"`; falling back to mgcv/gamlss."
    ))
  }

  if (!requireNamespace("scGLM", quietly = TRUE)) {
    return(unsupported("`use_scglm` requested but package 'scGLM' is not installed."))
  }
  if (!.scdesign3_is_intercept_only_formula(sigma_formula)) {
    return(unsupported("The scGLM backend currently supports only `sigma_formula = \"1\"`, except saturated categorical NB GAMLSS models."))
  }
  if (!family_key %in% c("binomial", "poisson", "gaussian", "nb")) {
    return(unsupported("The scGLM backend supports binomial, poisson, gaussian, and nb families."))
  }
  if (identical(family_key, "nb") && !requireNamespace("mgcv", quietly = TRUE)) {
    return(unsupported("The scGLM NB backend requires package 'mgcv'."))
  }

  fit_idx <- rep(TRUE, length(feature_names))
  if (!is.null(filtered_gene)) {
    fit_idx <- !feature_names %in% filtered_gene
  }
  if (!any(fit_idx)) {
    return(lapply(feature_names, function(gene) {
      .scdesign3_filtered_marginal(gene, trace = trace)
    }))
  }

  counts_fit <- count_mat[, fit_idx, drop = FALSE]
  if (!methods::is(counts_fit, "sparseMatrix")) {
    counts_fit <- as.matrix(counts_fit)
    storage.mode(counts_fit) <- "double"
  }
  rhs_formula <- stats::as.formula(paste0("~", mu_formula))
  full_formula <- stats::as.formula(paste0(predictor, "~", mu_formula))
  family_obj <- .scdesign3_scglm_family(family_key)
  method <- .scdesign3_scglm_method(
    method = scglm_method,
    family_key = family_key,
    formula = rhs_formula,
    data = dat_cov
  )
  if (identical(method, "unsupported")) {
    return(unsupported(paste0(
      "The automatic scGLM backend uses NB only for intercept-only or all-categorical formulas. ",
      "Set `use_scglm = \"always\"` and `scglm_method = \"newton_stein\"` to force an NB smooth fit."
    )))
  }
  scglm_fun <- function(name) {
    getExportedValue("scGLM", name)
  }

  start_time <- Sys.time()
  fit <- switch(
    method,
    categorical_closed_form = scglm_fun("scglm_categorical_closed_form_matrix")(
      counts = counts_fit,
      formula = rhs_formula,
      data = dat_cov,
      family = family_obj,
      batch_size = scglm_batch_size,
      n_cores = n_cores
    ),
    categorical_irls = scglm_fun("scglm_categorical_irls_matrix")(
      counts = counts_fit,
      formula = rhs_formula,
      data = dat_cov,
      family = family_obj,
      batch_size = scglm_batch_size
    ),
    newton_stein = scglm_fun("scglm_newton_stein_matrix")(
      counts = counts_fit,
      formula = rhs_formula,
      data = dat_cov,
      family = family_obj,
      batch_size = scglm_batch_size
    ),
    irls = scglm_fun("scglm_irls_matrix")(
      counts = counts_fit,
      formula = rhs_formula,
      data = dat_cov,
      family = family_obj,
      batch_size = scglm_batch_size
    ),
    stop("Unsupported scGLM method: ", method, call. = FALSE)
  )
  elapsed <- as.numeric(Sys.time() - start_time)

  x <- .scdesign3_scglm_training_x(fit)
  basis_object <- .scdesign3_scglm_basis_object(full_formula, dat_cov)
  rhs_terms <- stats::delete.response(stats::terms(rhs_formula, data = dat_cov))
  shared <- new.env(parent = emptyenv())
  shared$x <- x
  shared$data <- as.data.frame(dat_cov)
  shared$basis_object <- basis_object
  shared$design_names <- colnames(x)
  shared$rhs_terms <- rhs_terms
  shared$xlevels <- .scdesign3_xlevels(rhs_terms, dat_cov)
  shared$contrasts <- attr(x, "contrasts")
  cell_names <- rownames(count_mat)
  theta <- fit$theta
  if (is.null(theta)) {
    theta <- rep(NA_real_, ncol(counts_fit))
  }

  out <- vector("list", length(feature_names))
  names(out) <- feature_names
  fit_col <- 0L
  for (j in seq_along(feature_names)) {
    gene <- feature_names[j]
    if (!fit_idx[j]) {
      out[[gene]] <- .scdesign3_filtered_marginal(gene, trace = trace)
      next
    }

    fit_col <- fit_col + 1L
    y <- as.numeric(counts_fit[, fit_col])
    coef <- fit$coefficients[, fit_col]
    eta <- as.vector(x %*% coef)
    mu <- as.vector(fit$family$linkinv(eta))
    names(mu) <- cell_names
    names(eta) <- cell_names
    names(y) <- cell_names

    theta_g <- if (identical(family_key, "nb")) theta[fit_col] else NA_real_
    model <- .scdesign3_new_scglm_model(
      formula = full_formula,
      rhs_formula = rhs_formula,
      family_key = family_key,
      coefficients = coef,
      x = x,
      y = y,
      fitted = mu,
      eta = eta,
      theta = theta_g,
      data = dat_cov,
      shared = shared,
      backend = paste0("scGLM::", method),
      feature = gene
    )
    out[[gene]] <- if (isTRUE(trace)) {
      list(fit = model, warning = list(), time = c(elapsed, NA_real_), removed_cell = NA)
    } else {
      list(fit = model, removed_cell = NA)
    }
  }

  out
}

.scdesign3_fit_categorical_gamlss <- function(count_mat,
                                              dat_cov,
                                              filtered_gene,
                                              feature_names,
                                              family_key,
                                              mu_formula,
                                              sigma_formula,
                                              predictor,
                                              use_scglm,
                                              trace) {
  if (identical(family_key, "nb")) {
    return(.scdesign3_fit_categorical_nbi_gamlss(
      count_mat = count_mat,
      dat_cov = dat_cov,
      filtered_gene = filtered_gene,
      feature_names = feature_names,
      family_key = family_key,
      mu_formula = mu_formula,
      sigma_formula = sigma_formula,
      predictor = predictor,
      use_scglm = use_scglm,
      trace = trace
    ))
  }

  if (!family_key %in% c("binomial", "poisson", "gaussian", "zip", "zinb")) {
    return(NULL)
  }

  unsupported <- function(reason) {
    if (identical(use_scglm, "always")) {
      stop(reason, call. = FALSE)
    }
    NULL
  }

  sigma_intercept <- .scdesign3_is_intercept_only_formula(sigma_formula)
  use_fast <- switch(
    family_key,
    binomial = !sigma_intercept,
    poisson = !sigma_intercept,
    gaussian = !sigma_intercept,
    zip = TRUE,
    zinb = !sigma_intercept,
    FALSE
  )
  if (!isTRUE(use_fast)) {
    return(NULL)
  }

  mu_setup <- .scdesign3_saturated_categorical_setup(mu_formula, dat_cov)
  if (is.null(mu_setup)) {
    return(unsupported(
      "The fast categorical GAMLSS backend requires a saturated all-categorical mu formula."
    ))
  }

  if (family_key %in% c("gaussian", "zinb")) {
    sigma_setup <- .scdesign3_saturated_categorical_setup(sigma_formula, dat_cov)
    if (is.null(sigma_setup) ||
        !identical(mu_setup$design_names, sigma_setup$design_names) ||
        !identical(mu_setup$group, sigma_setup$group)) {
      return(unsupported(
        "The fast categorical GAMLSS backend currently requires identical saturated categorical mu and sigma formulas for this family."
      ))
    }
  }

  fit_idx <- rep(TRUE, length(feature_names))
  if (!is.null(filtered_gene)) {
    fit_idx <- !feature_names %in% filtered_gene
  }
  if (!any(fit_idx)) {
    return(lapply(feature_names, function(gene) {
      .scdesign3_filtered_marginal(gene, trace = trace)
    }))
  }

  start_time <- Sys.time()
  counts_fit <- count_mat[, fit_idx, drop = FALSE]
  group <- mu_setup$group
  params <- .scdesign3_categorical_family_params(
    counts_fit = counts_fit,
    group = group,
    n_groups = length(mu_setup$group_signatures),
    family_key = family_key
  )
  elapsed <- as.numeric(Sys.time() - start_time)

  shared <- new.env(parent = emptyenv())
  shared$rhs_formula <- mu_setup$rhs_formula
  shared$rhs_terms <- mu_setup$rhs_terms
  shared$xlevels <- mu_setup$xlevels
  shared$contrasts <- mu_setup$contrasts
  shared$design_names <- mu_setup$design_names
  shared$group_signatures <- mu_setup$group_signatures
  shared$train_group <- group
  shared$data <- as.data.frame(dat_cov)

  cell_names <- rownames(count_mat)
  out <- vector("list", length(feature_names))
  names(out) <- feature_names
  fit_col <- 0L
  for (j in seq_along(feature_names)) {
    gene <- feature_names[j]
    if (!fit_idx[j]) {
      out[[gene]] <- .scdesign3_filtered_marginal(gene, trace = trace)
      next
    }

    fit_col <- fit_col + 1L
    y <- as.numeric(counts_fit[, fit_col])
    names(y) <- cell_names
    model <- .scdesign3_new_categorical_dist_model(
      formula = stats::as.formula(paste0(predictor, "~", mu_formula)),
      sigma_formula = stats::as.formula(paste0("~", sigma_formula)),
      family_key = family_key,
      y = y,
      mu_by_group = params$mu[, fit_col],
      sigma_by_group = if (is.null(params$sigma)) NULL else params$sigma[, fit_col],
      nu_by_group = if (is.null(params$nu)) NULL else params$nu[, fit_col],
      df = params$df,
      shared = shared,
      feature = gene
    )
    out[[gene]] <- if (isTRUE(trace)) {
      list(
        fit = model,
        warning = list(),
        time = c(NA_real_, elapsed),
        removed_cell = NA
      )
    } else {
      list(fit = model, removed_cell = NA)
    }
  }

  out
}

.scdesign3_fit_shared_fixed_smooth <- function(count_mat,
                                               dat_cov,
                                               filtered_gene,
                                               feature_names,
                                               family_key,
                                               mu_formula,
                                               sigma_formula,
                                               predictor,
                                               use_scglm,
                                               scglm_method,
                                               scglm_batch_size,
                                               n_cores,
                                               trace) {
  if (!.scdesign3_formula_has_smooth(mu_formula)) {
    return(NULL)
  }

  unsupported <- function(reason) {
    if (identical(use_scglm, "always")) {
      stop(reason, call. = FALSE)
    }
    NULL
  }

  if (!.scdesign3_is_intercept_only_formula(sigma_formula)) {
    return(unsupported("The shared fixed-smooth backend currently requires `sigma_formula = \"1\"`."))
  }
  if (!family_key %in% c("poisson", "binomial", "gaussian")) {
    return(unsupported("The shared fixed-smooth backend currently supports poisson, binomial, and gaussian families."))
  }
  if (!scglm_method %in% c("auto", "irls")) {
    return(unsupported("The shared fixed-smooth backend currently supports `scglm_method = \"auto\"` or `\"irls\"`."))
  }

  setup <- .scdesign3_fixed_smooth_setup(
    mu_formula = mu_formula,
    predictor = predictor,
    dat_cov = dat_cov,
    family_key = family_key
  )
  if (is.null(setup)) {
    return(unsupported(
      "The smooth formula is not fixed. Use `fx = TRUE` or non-negative `sp = ...` inside the smooth term to use the shared fixed-smooth backend."
    ))
  }

  fit_idx <- rep(TRUE, length(feature_names))
  if (!is.null(filtered_gene)) {
    fit_idx <- !feature_names %in% filtered_gene
  }
  if (!any(fit_idx)) {
    return(lapply(feature_names, function(gene) {
      .scdesign3_filtered_marginal(gene, trace = trace)
    }))
  }

  counts_fit <- as.matrix(count_mat[, fit_idx, drop = FALSE])
  storage.mode(counts_fit) <- "double"
  family_obj <- .scdesign3_fixed_smooth_family(family_key)

  start_time <- Sys.time()
  fit <- if (!isTRUE(setup$has_penalty) &&
             identical(family_key, "poisson") &&
             requireNamespace("scGLM", quietly = TRUE)) {
    .scdesign3_shared_unpenalized_poisson_irls(
      x = setup$x,
      y = counts_fit,
      maxit = 100L,
      tol = 1e-10,
      batch_size = scglm_batch_size,
      n_cores = n_cores
    )
  } else {
    .scdesign3_shared_penalized_irls(
      x = setup$x,
      y = counts_fit,
      family = family_obj,
      penalty = setup$penalty,
      maxit = 100L,
      tol = 1e-10,
      batch_size = scglm_batch_size,
      n_cores = n_cores
    )
  }
  elapsed <- as.numeric(Sys.time() - start_time)

  rhs_formula <- stats::as.formula(paste0("~", mu_formula))
  full_formula <- stats::as.formula(paste0(predictor, "~", mu_formula))
  basis_object <- .scdesign3_scglm_basis_object(full_formula, dat_cov)
  rhs_terms <- stats::delete.response(stats::terms(rhs_formula, data = dat_cov))
  shared <- new.env(parent = emptyenv())
  shared$x <- setup$x
  shared$data <- as.data.frame(dat_cov)
  shared$basis_object <- basis_object
  shared$design_names <- colnames(setup$x)
  shared$rhs_terms <- rhs_terms
  shared$xlevels <- .scdesign3_xlevels(rhs_terms, dat_cov)
  shared$contrasts <- NULL
  shared$penalty <- setup$penalty

  cell_names <- rownames(count_mat)
  out <- vector("list", length(feature_names))
  names(out) <- feature_names
  fit_col <- 0L
  for (j in seq_along(feature_names)) {
    gene <- feature_names[j]
    if (!fit_idx[j]) {
      out[[gene]] <- .scdesign3_filtered_marginal(gene, trace = trace)
      next
    }

    fit_col <- fit_col + 1L
    y <- counts_fit[, fit_col]
    coef <- fit$coefficients[, fit_col]
    eta <- fit$eta[, fit_col]
    mu <- fit$mu[, fit_col]
    names(mu) <- cell_names
    names(eta) <- cell_names
    names(y) <- cell_names

    model <- .scdesign3_new_scglm_model(
      formula = full_formula,
      rhs_formula = rhs_formula,
      family_key = family_key,
      coefficients = coef,
      x = setup$x,
      y = y,
      fitted = mu,
      eta = eta,
      theta = NA_real_,
      data = dat_cov,
      shared = shared,
      backend = setup$backend,
      feature = gene
    )
    model$df <- fit$edf[fit_col] + as.integer(identical(family_key, "gaussian"))
    model$edf <- fit$edf[fit_col]
    model$edf1 <- fit$edf[fit_col]
    model$edf2 <- fit$edf[fit_col]
    model$df.residual <- max(length(y) - fit$edf[fit_col], 0)
    model$logLik <- fit$logLik[fit_col]
    if (identical(family_key, "gaussian")) {
      model$sig2 <- fit$sig2[fit_col]
    }

    out[[gene]] <- if (isTRUE(trace)) {
      list(fit = model, warning = list(), time = c(elapsed, NA_real_), removed_cell = NA)
    } else {
      list(fit = model, removed_cell = NA)
    }
  }

  out
}

.scdesign3_fixed_smooth_setup <- function(mu_formula,
                                          predictor,
                                          dat_cov,
                                          family_key) {
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(NULL)
  }
  full_formula <- stats::as.formula(paste0(predictor, "~", mu_formula))
  dat_template <- as.data.frame(dat_cov)
  dat_template[[predictor]] <- .scdesign3_mgcv_placeholder_y(nrow(dat_template), family_key)
  G <- tryCatch(
    mgcv::gam(
      formula = full_formula,
      data = dat_template,
      family = .scdesign3_fixed_smooth_family(family_key),
      fit = FALSE,
      method = "REML"
    ),
    error = function(e) NULL
  )
  if (is.null(G) || is.null(G$X) || length(G$smooth) == 0L) {
    return(NULL)
  }

  penalty <- matrix(0, nrow = ncol(G$X), ncol = ncol(G$X))
  has_penalty <- FALSE
  for (sm in G$smooth) {
    if (length(sm$S) == 0L) {
      next
    }
    sp <- sm$sp
    if (is.null(sp) || length(sp) != length(sm$S) ||
        any(!is.finite(sp)) || any(sp < 0)) {
      return(NULL)
    }
    for (i in seq_along(sm$S)) {
      s_i <- sm$S[[i]]
      idx <- seq.int(sm$first.para, length.out = nrow(s_i))
      penalty[idx, idx] <- penalty[idx, idx, drop = FALSE] + sp[i] * s_i
      has_penalty <- TRUE
    }
  }

  x <- G$X
  colnames(x) <- .scdesign3_mgcv_design_colnames(G)
  list(
    x = x,
    penalty = penalty,
    has_penalty = has_penalty,
    backend = if (has_penalty) {
      "scDesign3::shared_fixed_penalty_smooth_irls"
    } else {
      "scDesign3::shared_fixed_basis_smooth_irls"
    }
  )
}

.scdesign3_fixed_smooth_family <- function(family_key) {
  switch(
    family_key,
    poisson = stats::poisson(),
    binomial = stats::binomial(),
    gaussian = stats::gaussian(),
    stop("Unsupported fixed-smooth family: ", family_key, call. = FALSE)
  )
}

.scdesign3_mgcv_design_colnames <- function(gam_obj) {
  x_names <- colnames(gam_obj$X)
  if (is.null(x_names)) {
    x_names <- rep.int("", ncol(gam_obj$X))
  }
  x_names[is.na(x_names)] <- ""
  for (sm in gam_obj$smooth) {
    idx <- seq.int(sm$first.para, sm$last.para)
    blank <- !nzchar(x_names[idx])
    x_names[idx][blank] <- paste0(sm$label, ".", seq_len(length(idx))[blank])
  }
  blanks <- which(!nzchar(x_names))
  if (length(blanks) > 0L) {
    x_names[blanks] <- paste0("basis.", blanks)
  }
  make.unique(x_names)
}

.scdesign3_shared_unpenalized_poisson_irls <- function(x,
                                                       y,
                                                       maxit = 100L,
                                                       tol = 1e-10,
                                                       batch_size = 256L,
                                                       n_cores = 1L) {
  beta_init <- .scdesign3_initial_beta(x, y, stats::poisson())
  blocks <- .scdesign3_gene_blocks(ncol(y), batch_size = batch_size)
  fit_block <- function(idx) {
    fit <- getFromNamespace("scglm_poisson_irls_cpp", "scGLM")(
      x = x,
      y = y[, idx, drop = FALSE],
      beta_init = beta_init[, idx, drop = FALSE],
      maxit = as.integer(maxit),
      tol = tol
    )
    list(
      idx = idx,
      coefficients = fit$coefficients,
      iterations = as.integer(fit$iterations),
      converged = as.logical(fit$converged)
    )
  }
  block_fit <- .scdesign3_parallel_lapply(blocks, fit_block, n_cores = n_cores)
  beta <- matrix(0, nrow = ncol(x), ncol = ncol(y))
  iterations <- integer(ncol(y))
  converged <- logical(ncol(y))
  for (res in block_fit) {
    beta[, res$idx] <- res$coefficients
    iterations[res$idx] <- res$iterations
    converged[res$idx] <- res$converged
  }
  rownames(beta) <- colnames(x)
  colnames(beta) <- colnames(y)
  eta <- pmax(pmin(x %*% beta, 30), -30)
  mu <- exp(eta)
  ll <- colSums(stats::dpois(y, lambda = pmax(mu, .Machine$double.eps), log = TRUE))
  list(
    coefficients = beta,
    eta = eta,
    mu = mu,
    edf = rep(qr(x)$rank, ncol(y)),
    sig2 = rep(1, ncol(y)),
    logLik = ll,
    iterations = iterations,
    converged = converged
  )
}

.scdesign3_shared_penalized_irls <- function(x,
                                             y,
                                             family,
                                             penalty,
                                             maxit = 100L,
                                             tol = 1e-10,
                                             batch_size = 256L,
                                             n_cores = 1L) {
  n <- nrow(x)
  p <- ncol(x)
  g <- ncol(y)
  penalty <- as.matrix(penalty)
  if (identical(family$family, "gaussian") && identical(family$link, "identity")) {
    h <- crossprod(x) + penalty
    beta <- .scdesign3_safe_solve(h, crossprod(x, y))
    eta <- x %*% beta
    mu <- eta
    edf <- rep(.scdesign3_hat_trace(h, crossprod(x)), g)
    rss <- colSums((y - mu)^2)
    sig2 <- rss / pmax(n - edf, 1e-8)
    sig2_mle <- rss / n
    ll <- vapply(seq_len(g), function(j) {
      sum(stats::dnorm(y[, j], mean = mu[, j], sd = sqrt(pmax(sig2_mle[j], .Machine$double.eps)), log = TRUE))
    }, numeric(1))
    rownames(beta) <- colnames(x)
    colnames(beta) <- colnames(y)
    return(list(
      coefficients = beta,
      eta = eta,
      mu = mu,
      edf = edf,
      sig2 = sig2,
      logLik = ll,
      iterations = 1L,
      converged = rep(TRUE, g)
    ))
  }

  beta <- .scdesign3_initial_beta(x, y, family)
  converged <- rep(FALSE, g)
  for (iter in seq_len(maxit)) {
    beta_old <- beta
    eta <- pmax(pmin(x %*% beta, 30), -30)
    mu <- family$linkinv(eta)
    mu_eta <- pmax(abs(family$mu.eta(eta)), .Machine$double.eps)
    variance <- pmax(family$variance(mu), .Machine$double.eps)
    w <- pmax((mu_eta^2) / variance, .Machine$double.eps)
    z <- eta + (y - mu) / mu_eta

    blocks <- .scdesign3_gene_blocks(g, batch_size = batch_size)
    fit_block <- function(block) {
      beta_block <- matrix(0, nrow = p, ncol = length(block))
      for (k in seq_along(block)) {
        j <- block[k]
        xw <- x * sqrt(w[, j])
        h <- crossprod(xw) + penalty
        rhs <- crossprod(x, w[, j] * z[, j])
        beta_block[, k] <- .scdesign3_safe_solve(h, rhs)
      }
      list(idx = block, beta = beta_block)
    }
    block_fit <- .scdesign3_parallel_lapply(blocks, fit_block, n_cores = n_cores)
    for (res in block_fit) {
      beta[, res$idx] <- res$beta
    }
    delta <- max(abs(beta - beta_old) / pmax(abs(beta_old), 0.1))
    if (is.finite(delta) && delta < tol) {
      converged[] <- TRUE
      break
    }
  }

  eta <- pmax(pmin(x %*% beta, 30), -30)
  mu <- family$linkinv(eta)
  mu_eta <- pmax(abs(family$mu.eta(eta)), .Machine$double.eps)
  variance <- pmax(family$variance(mu), .Machine$double.eps)
  w <- pmax((mu_eta^2) / variance, .Machine$double.eps)
  edf <- numeric(g)
  ll <- numeric(g)
  blocks <- .scdesign3_gene_blocks(g, batch_size = batch_size)
  summary_block <- function(block) {
    edf_block <- numeric(length(block))
    ll_block <- numeric(length(block))
    for (k in seq_along(block)) {
      j <- block[k]
      xw <- x * sqrt(w[, j])
      xtwx <- crossprod(xw)
      h <- xtwx + penalty
      edf_block[k] <- .scdesign3_hat_trace(h, xtwx)
      ll_block[k] <- switch(
        family$family,
        poisson = sum(stats::dpois(y[, j], lambda = pmax(mu[, j], .Machine$double.eps), log = TRUE)),
        binomial = sum(stats::dbinom(y[, j], size = 1, prob = pmin(pmax(mu[, j], .Machine$double.eps), 1 - .Machine$double.eps), log = TRUE)),
        sum(family$dev.resids(y[, j], mu[, j], rep.int(1, n))) / -2
      )
    }
    list(idx = block, edf = edf_block, logLik = ll_block)
  }
  block_summary <- .scdesign3_parallel_lapply(blocks, summary_block, n_cores = n_cores)
  for (res in block_summary) {
    edf[res$idx] <- res$edf
    ll[res$idx] <- res$logLik
  }

  rownames(beta) <- colnames(x)
  colnames(beta) <- colnames(y)
  list(
    coefficients = beta,
    eta = eta,
    mu = mu,
    edf = edf,
    sig2 = rep(1, g),
    logLik = ll,
    iterations = iter,
    converged = converged
  )
}

.scdesign3_sanitize_n_cores <- function(n_cores, n_tasks = Inf) {
  n_cores <- suppressWarnings(as.integer(n_cores[1L]))
  if (!is.finite(n_cores) || n_cores < 1L) {
    n_cores <- 1L
  }
  max(1L, min(n_cores, as.integer(n_tasks)))
}

.scdesign3_gene_blocks <- function(n_genes, batch_size = 256L) {
  if (n_genes <= 0L) {
    return(list())
  }
  batch_size <- suppressWarnings(as.integer(batch_size[1L]))
  if (!is.finite(batch_size) || batch_size < 1L) {
    batch_size <- 256L
  }
  split(seq_len(n_genes), ceiling(seq_len(n_genes) / batch_size))
}

.scdesign3_parallel_lapply <- function(x, fun, n_cores = 1L) {
  n_cores <- .scdesign3_sanitize_n_cores(n_cores, length(x))
  if (n_cores <= 1L || length(x) <= 1L || .Platform$OS.type == "windows") {
    return(lapply(x, fun))
  }
  parallel::mclapply(x, fun, mc.cores = n_cores, mc.preschedule = TRUE)
}

.scdesign3_initial_beta <- function(x, y, family) {
  g <- ncol(y)
  beta <- matrix(0, nrow = ncol(x), ncol = g)
  rownames(beta) <- colnames(x)
  intercept <- match("(Intercept)", colnames(x), nomatch = 1L)
  ybar <- colMeans(y)
  beta[intercept, ] <- switch(
    family$family,
    poisson = log(pmax(ybar, 0.1)),
    binomial = stats::qlogis(pmin(pmax(ybar, 1e-4), 1 - 1e-4)),
    gaussian = ybar,
    ybar
  )
  beta
}

.scdesign3_safe_solve <- function(a, b) {
  out <- tryCatch(
    solve(a, b),
    error = function(e) qr.solve(a, b)
  )
  out[!is.finite(out)] <- 0
  out
}

.scdesign3_hat_trace <- function(h, xtwx) {
  sum(diag(.scdesign3_safe_solve(h, xtwx)))
}

.scdesign3_fit_categorical_nbi_gamlss <- function(count_mat,
                                                  dat_cov,
                                                  filtered_gene,
                                                  feature_names,
                                                  family_key,
                                                  mu_formula,
                                                  sigma_formula,
                                                  predictor,
                                                  use_scglm,
                                                  trace) {
  if (!identical(family_key, "nb") ||
      .scdesign3_is_intercept_only_formula(sigma_formula)) {
    return(NULL)
  }

  unsupported <- function(reason) {
    if (identical(use_scglm, "always")) {
      stop(reason, call. = FALSE)
    }
    NULL
  }

  mu_setup <- .scdesign3_saturated_categorical_setup(mu_formula, dat_cov)
  sigma_setup <- .scdesign3_saturated_categorical_setup(sigma_formula, dat_cov)
  if (is.null(mu_setup) || is.null(sigma_setup)) {
    return(unsupported(
      "The fast NB GAMLSS backend requires saturated all-categorical mu and sigma formulas."
    ))
  }
  if (!identical(mu_setup$design_names, sigma_setup$design_names) ||
      !identical(mu_setup$group, sigma_setup$group)) {
    return(unsupported(
      "The fast NB GAMLSS backend currently requires identical saturated categorical mu and sigma formulas."
    ))
  }

  fit_idx <- rep(TRUE, length(feature_names))
  if (!is.null(filtered_gene)) {
    fit_idx <- !feature_names %in% filtered_gene
  }
  if (!any(fit_idx)) {
    return(lapply(feature_names, function(gene) {
      .scdesign3_filtered_marginal(gene, trace = trace)
    }))
  }

  start_time <- Sys.time()
  counts_fit <- count_mat[, fit_idx, drop = FALSE]
  group <- mu_setup$group
  n_by_group <- as.numeric(tabulate(group, nbins = length(mu_setup$group_signatures)))
  z <- Matrix::sparseMatrix(
    i = seq_along(group),
    j = group,
    x = 1,
    dims = c(length(group), length(mu_setup$group_signatures))
  )

  sum_y <- as.matrix(Matrix::crossprod(z, counts_fit))
  counts_sq <- counts_fit^2
  sum_y2 <- as.matrix(Matrix::crossprod(z, counts_sq))
  mu_by_group <- sweep(sum_y, 1L, n_by_group, "/")
  second_moment <- sweep(sum_y2, 1L, n_by_group, "/")
  var_by_group <- pmax(second_moment - mu_by_group^2, 0)
  if (any(n_by_group > 1)) {
    var_by_group <- sweep(var_by_group, 1L, pmax(n_by_group - 1, 1) / pmax(n_by_group, 1), "/")
  }
  sigma_start <- pmax((var_by_group - mu_by_group) / pmax(mu_by_group^2, .Machine$double.eps), 1e-8)
  sigma_start[!is.finite(sigma_start)] <- 1e-8
  sigma_by_group <- .scdesign3_nbi_sigma_group_mle(
    counts_fit = counts_fit,
    group = group,
    mu_by_group = mu_by_group,
    sigma_start = sigma_start
  )
  sigma_constant <- .scdesign3_nbi_sigma_constant_mle(
    counts_fit = counts_fit,
    group = group,
    mu_by_group = mu_by_group
  )
  zero_mu <- mu_by_group <= 0 | !is.finite(mu_by_group)
  if (any(zero_mu)) {
    sigma_by_group[zero_mu] <- sigma_constant[col(zero_mu)[zero_mu]]
  }
  bad_sigma <- apply(
    !is.finite(sigma_by_group) | sigma_by_group > 1000,
    2L,
    any
  )
  if (any(bad_sigma)) {
    sigma_by_group[, bad_sigma] <- matrix(
      rep(sigma_constant[bad_sigma], each = nrow(sigma_by_group)),
      nrow = nrow(sigma_by_group)
    )
  }
  elapsed <- as.numeric(Sys.time() - start_time)

  shared <- new.env(parent = emptyenv())
  shared$rhs_formula <- mu_setup$rhs_formula
  shared$rhs_terms <- mu_setup$rhs_terms
  shared$xlevels <- mu_setup$xlevels
  shared$contrasts <- mu_setup$contrasts
  shared$design_names <- mu_setup$design_names
  shared$group_signatures <- mu_setup$group_signatures
  shared$train_group <- group
  shared$data <- as.data.frame(dat_cov)

  cell_names <- rownames(count_mat)
  out <- vector("list", length(feature_names))
  names(out) <- feature_names
  fit_col <- 0L
  for (j in seq_along(feature_names)) {
    gene <- feature_names[j]
    if (!fit_idx[j]) {
      out[[gene]] <- .scdesign3_filtered_marginal(gene, trace = trace)
      next
    }

    fit_col <- fit_col + 1L
    y <- as.numeric(counts_fit[, fit_col])
    names(y) <- cell_names
    model <- .scdesign3_new_categorical_nbi_model(
      formula = stats::as.formula(paste0(predictor, "~", mu_formula)),
      sigma_formula = stats::as.formula(paste0("~", sigma_formula)),
      y = y,
      mu_by_group = mu_by_group[, fit_col],
      sigma_by_group = sigma_by_group[, fit_col],
      shared = shared,
      feature = gene
    )
    out[[gene]] <- if (isTRUE(trace)) {
      list(
        fit = model,
        warning = list(),
        time = c(NA_real_, elapsed),
        removed_cell = NA
      )
    } else {
      list(fit = model, removed_cell = NA)
    }
  }

  out
}

.scdesign3_saturated_categorical_setup <- function(formula_string, data) {
  rhs_formula <- stats::as.formula(paste0("~", formula_string))
  if (!.scdesign3_formula_is_all_categorical(rhs_formula, data)) {
    return(NULL)
  }

  rhs_terms <- stats::terms(rhs_formula, data = data)
  rhs_formula <- .scdesign3_compact_formula(rhs_formula)
  rhs_terms <- .scdesign3_compact_terms(rhs_terms)
  rhs_frame <- stats::model.frame(rhs_terms, data = data, na.action = stats::na.pass)
  if (anyNA(rhs_frame)) {
    return(NULL)
  }
  x <- stats::model.matrix(rhs_terms, data = rhs_frame)
  design_signatures <- .scdesign3_model_matrix_signature(x)
  group_signatures <- unique(design_signatures)
  group <- match(design_signatures, group_signatures)

  rank_x <- qr(x)$rank
  if (!identical(rank_x, length(group_signatures))) {
    return(NULL)
  }

  list(
    rhs_formula = rhs_formula,
    rhs_terms = rhs_terms,
    xlevels = .scdesign3_xlevels(rhs_terms, data),
    contrasts = attr(x, "contrasts"),
    design_names = colnames(x),
    group = group,
    group_signatures = group_signatures
  )
}

.scdesign3_compact_formula <- function(formula) {
  environment(formula) <- baseenv()
  formula
}

.scdesign3_compact_terms <- function(terms) {
  attr(terms, ".Environment") <- baseenv()
  terms
}

.scdesign3_model_matrix_signature <- function(x) {
  if (is.null(dim(x)) || nrow(x) == 0L) {
    return(character())
  }
  apply(
    x,
    1L,
    function(row) paste(format(signif(row, 14L), scientific = FALSE, trim = TRUE), collapse = "\r")
  )
}

.scdesign3_nbi_sigma_group_mle <- function(counts_fit,
                                           group,
                                           mu_by_group,
                                           sigma_start,
                                           lower = 1e-8,
                                           upper = 1e4) {
  out <- sigma_start
  log_lower <- log(lower)
  log_upper <- log(upper)
  group_index <- split(seq_along(group), group)
  for (j in seq_len(ncol(counts_fit))) {
    y_col <- as.numeric(counts_fit[, j])
    for (k in seq_len(nrow(mu_by_group))) {
      idx <- group_index[[as.character(k)]]
      mu <- mu_by_group[k, j]
      if (!is.finite(mu) || mu <= 0 || length(idx) < 2L) {
        out[k, j] <- lower
        next
      }
      y <- y_col[idx]
      if (stats::var(y) <= mu) {
        out[k, j] <- lower
        next
      }
      tab <- table(y)
      yy <- as.numeric(names(tab))
      ww <- as.numeric(tab)
      objective <- function(log_sigma) {
        sigma <- exp(log_sigma)
        -sum(ww * stats::dnbinom(yy, size = 1 / sigma, mu = mu, log = TRUE))
      }
      fit <- tryCatch(
        stats::optimize(objective, interval = c(log_lower, log_upper)),
        error = function(e) NULL
      )
      if (is.null(fit) || !is.finite(fit$minimum)) {
        out[k, j] <- min(max(sigma_start[k, j], lower), upper)
      } else {
        out[k, j] <- min(max(exp(fit$minimum), lower), upper)
      }
    }
  }
  rownames(out) <- rownames(mu_by_group)
  colnames(out) <- colnames(mu_by_group)
  out
}

.scdesign3_nbi_sigma_constant_mle <- function(counts_fit,
                                              group,
                                              mu_by_group,
                                              lower = 1e-8,
                                              upper = 1e4) {
  out <- rep(lower, ncol(counts_fit))
  log_lower <- log(lower)
  log_upper <- log(upper)
  for (j in seq_len(ncol(counts_fit))) {
    y <- as.numeric(counts_fit[, j])
    mu <- pmax(mu_by_group[group, j], .Machine$double.eps)
    if (stats::var(y) <= mean(mu)) {
      out[j] <- lower
      next
    }
    objective <- function(log_sigma) {
      sigma <- exp(log_sigma)
      -sum(stats::dnbinom(y, size = 1 / sigma, mu = mu, log = TRUE))
    }
    fit <- tryCatch(
      stats::optimize(objective, interval = c(log_lower, log_upper)),
      error = function(e) NULL
    )
    if (!is.null(fit) && is.finite(fit$minimum)) {
      out[j] <- min(max(exp(fit$minimum), lower), upper)
    }
  }
  out
}

.scdesign3_categorical_family_params <- function(counts_fit,
                                                 group,
                                                 n_groups,
                                                 family_key) {
  z <- Matrix::sparseMatrix(
    i = seq_along(group),
    j = group,
    x = 1,
    dims = c(length(group), n_groups)
  )
  n_by_group <- as.numeric(tabulate(group, nbins = n_groups))
  sum_y <- as.matrix(Matrix::crossprod(z, counts_fit))
  mu_obs <- sweep(sum_y, 1L, n_by_group, "/")

  if (family_key == "poisson") {
    return(list(mu = mu_obs, sigma = NULL, nu = NULL, df = n_groups))
  }

  if (family_key == "binomial") {
    mu <- pmin(pmax(mu_obs, 0), 1)
    return(list(mu = mu, sigma = NULL, nu = NULL, df = n_groups))
  }

  if (family_key == "gaussian") {
    sum_y2 <- as.matrix(Matrix::crossprod(z, counts_fit^2))
    second_moment <- sweep(sum_y2, 1L, n_by_group, "/")
    sigma <- sqrt(pmax(second_moment - mu_obs^2, .Machine$double.eps))
    return(list(mu = mu_obs, sigma = sigma, nu = NULL, df = 2L * n_groups))
  }

  if (family_key == "zip") {
    return(.scdesign3_categorical_zip_params(counts_fit, group, n_groups))
  }

  if (family_key == "zinb") {
    return(.scdesign3_categorical_zinb_params(counts_fit, group, n_groups))
  }

  stop("Unsupported categorical family: ", family_key, call. = FALSE)
}

.scdesign3_categorical_zip_params <- function(counts_fit, group, n_groups) {
  mu <- matrix(0, nrow = n_groups, ncol = ncol(counts_fit))
  sigma <- matrix(0, nrow = n_groups, ncol = ncol(counts_fit))
  group_index <- split(seq_along(group), group)
  for (j in seq_len(ncol(counts_fit))) {
    y_col <- as.numeric(counts_fit[, j])
    for (k in seq_len(n_groups)) {
      pars <- .scdesign3_zip_group_mle(y_col[group_index[[as.character(k)]]])
      mu[k, j] <- pars[["mu"]]
      sigma[k, j] <- pars[["sigma"]]
    }
  }
  rownames(mu) <- rownames(sigma) <- seq_len(n_groups)
  colnames(mu) <- colnames(sigma) <- colnames(counts_fit)
  list(mu = mu, sigma = sigma, nu = NULL, df = 2L * n_groups)
}

.scdesign3_categorical_zinb_params <- function(counts_fit, group, n_groups) {
  mu <- sigma <- nu <- matrix(0, nrow = n_groups, ncol = ncol(counts_fit))
  group_index <- split(seq_along(group), group)
  for (j in seq_len(ncol(counts_fit))) {
    y_col <- as.numeric(counts_fit[, j])
    for (k in seq_len(n_groups)) {
      pars <- .scdesign3_zinb_group_mle(y_col[group_index[[as.character(k)]]])
      mu[k, j] <- pars[["mu"]]
      sigma[k, j] <- pars[["sigma"]]
      nu[k, j] <- pars[["nu"]]
    }
  }
  rownames(mu) <- rownames(sigma) <- rownames(nu) <- seq_len(n_groups)
  colnames(mu) <- colnames(sigma) <- colnames(nu) <- colnames(counts_fit)
  list(mu = mu, sigma = sigma, nu = nu, df = 3L * n_groups)
}

.scdesign3_zip_group_mle <- function(y, eps = 1e-8) {
  y <- as.numeric(y)
  if (length(y) == 0L || all(!is.finite(y))) {
    return(c(mu = eps, sigma = eps))
  }
  y <- y[is.finite(y)]
  if (all(y <= 0)) {
    return(c(mu = eps, sigma = 1 - eps))
  }
  if (!any(y == 0)) {
    return(c(mu = max(mean(y), eps), sigma = eps))
  }

  tab <- table(y)
  yy <- as.numeric(names(tab))
  ww <- as.numeric(tab)
  mean_y <- mean(y)
  mu_start <- max(mean(y[y > 0]), mean_y, 0.1)
  pois_zero <- exp(-mu_start)
  sigma_start <- (mean(y == 0) - pois_zero) / max(1 - pois_zero, eps)
  sigma_start <- pmin(pmax(sigma_start, eps), 1 - eps)
  upper_mu <- max(max(y) * 20, mean_y * 50, 10) + 1

  objective <- function(par) {
    mu <- exp(par[1L])
    sigma <- stats::plogis(par[2L])
    -sum(ww * gamlss.dist::dZIP(yy, mu = mu, sigma = sigma, log = TRUE))
  }
  fit <- tryCatch(
    suppressWarnings(stats::optim(
      par = c(log(mu_start), stats::qlogis(sigma_start)),
      fn = objective,
      method = "L-BFGS-B",
      lower = c(log(eps), -30),
      upper = c(log(upper_mu), 30)
    )),
    error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(fit$value)) {
    return(c(mu = mu_start, sigma = sigma_start))
  }
  c(mu = exp(fit$par[1L]), sigma = stats::plogis(fit$par[2L]))
}

.scdesign3_zinb_group_mle <- function(y, eps = 1e-8) {
  y <- as.numeric(y)
  if (length(y) == 0L || all(!is.finite(y))) {
    return(c(mu = eps, sigma = eps, nu = eps))
  }
  y <- y[is.finite(y)]
  if (all(y <= 0)) {
    return(c(mu = eps, sigma = eps, nu = 1 - eps))
  }

  tab <- table(y)
  yy <- as.numeric(names(tab))
  ww <- as.numeric(tab)
  mean_y <- mean(y)
  var_y <- stats::var(y)
  mu_start <- max(mean(y[y > 0]), mean_y, 0.1)
  sigma_start <- pmax((var_y - mean_y) / pmax(mean_y^2, eps), 1e-4)
  sigma_start <- pmin(pmax(sigma_start, eps), 1e4)
  nb_zero <- stats::dnbinom(0, size = 1 / sigma_start, mu = mu_start)
  nu_start <- (mean(y == 0) - nb_zero) / max(1 - nb_zero, eps)
  nu_start <- pmin(pmax(nu_start, eps), 1 - eps)
  upper_mu <- max(max(y) * 20, mean_y * 50, 10) + 1

  objective <- function(par) {
    mu <- exp(par[1L])
    sigma <- exp(par[2L])
    nu <- stats::plogis(par[3L])
    -sum(ww * gamlss.dist::dZINBI(yy, mu = mu, sigma = sigma, nu = nu, log = TRUE))
  }
  fit <- tryCatch(
    suppressWarnings(stats::optim(
      par = c(log(mu_start), log(sigma_start), stats::qlogis(nu_start)),
      fn = objective,
      method = "L-BFGS-B",
      lower = c(log(eps), log(eps), -30),
      upper = c(log(upper_mu), log(1e4), 30)
    )),
    error = function(e) NULL
  )
  if (is.null(fit) || !is.finite(fit$value)) {
    return(c(mu = mu_start, sigma = sigma_start, nu = nu_start))
  }
  c(
    mu = exp(fit$par[1L]),
    sigma = exp(fit$par[2L]),
    nu = stats::plogis(fit$par[3L])
  )
}

.scdesign3_new_categorical_dist_model <- function(formula,
                                                  sigma_formula,
                                                  family_key,
                                                  y,
                                                  mu_by_group,
                                                  sigma_by_group,
                                                  nu_by_group,
                                                  df,
                                                  shared,
                                                  feature) {
  formula <- .scdesign3_compact_formula(formula)
  sigma_formula <- .scdesign3_compact_formula(sigma_formula)
  names(mu_by_group) <- shared$group_signatures
  if (!is.null(sigma_by_group)) {
    names(sigma_by_group) <- shared$group_signatures
  }
  if (!is.null(nu_by_group)) {
    names(nu_by_group) <- shared$group_signatures
  }
  fitted <- mu_by_group[shared$group_signatures[shared$train_group]]
  sigma_fitted <- if (is.null(sigma_by_group)) {
    rep(NA_real_, length(fitted))
  } else {
    sigma_by_group[shared$group_signatures[shared$train_group]]
  }
  nu_fitted <- if (is.null(nu_by_group)) {
    rep(0, length(fitted))
  } else {
    nu_by_group[shared$group_signatures[shared$train_group]]
  }
  names(fitted) <- names(sigma_fitted) <- names(nu_fitted) <- names(y)
  ll <- .scdesign3_categorical_loglik(
    y = y,
    family_key = family_key,
    mu = fitted,
    sigma = sigma_fitted,
    nu = nu_fitted
  )

  structure(
    list(
      formula = formula,
      terms = .scdesign3_compact_terms(stats::terms(formula, data = shared$data)),
      mu.formula = formula,
      sigma.formula = sigma_formula,
      family = .scdesign3_categorical_family(family_key),
      family_key = family_key,
      y = y,
      mu_by_group = mu_by_group,
      sigma_by_group = sigma_by_group,
      nu_by_group = nu_by_group,
      logLik = ll,
      df = df,
      df.residual = max(length(y) - df, 0L),
      nobs = length(y),
      shared = shared,
      backend = paste0("scDesign3::categorical_", family_key, "_group_mle"),
      feature = feature
    ),
    class = "scdesign3_categorical_dist"
  )
}

.scdesign3_categorical_loglik <- function(y, family_key, mu, sigma = NULL, nu = NULL) {
  mu <- as.numeric(mu)
  sigma <- as.numeric(sigma)
  nu <- as.numeric(nu)
  switch(
    family_key,
    binomial = sum(stats::dbinom(y, size = 1, prob = pmin(pmax(mu, 0), 1), log = TRUE)),
    poisson = sum(stats::dpois(y, lambda = pmax(mu, 0), log = TRUE)),
    gaussian = sum(stats::dnorm(y, mean = mu, sd = pmax(sigma, .Machine$double.eps), log = TRUE)),
    zip = sum(gamlss.dist::dZIP(y, mu = pmax(mu, .Machine$double.eps), sigma = pmin(pmax(sigma, .Machine$double.eps), 1 - .Machine$double.eps), log = TRUE)),
    zinb = sum(gamlss.dist::dZINBI(y, mu = pmax(mu, .Machine$double.eps), sigma = pmax(sigma, .Machine$double.eps), nu = pmin(pmax(nu, .Machine$double.eps), 1 - .Machine$double.eps), log = TRUE)),
    stop("Unsupported categorical family log-likelihood: ", family_key, call. = FALSE)
  )
}

.scdesign3_categorical_family <- function(family_key) {
  list(
    family = switch(
      family_key,
      binomial = "BI",
      poisson = "PO",
      gaussian = "NO",
      zip = "ZIP",
      zinb = "ZINBI",
      family_key
    ),
    link = switch(
      family_key,
      binomial = "logit",
      gaussian = "identity",
      "log"
    )
  )
}

.scdesign3_new_categorical_nbi_model <- function(formula,
                                                 sigma_formula,
                                                 y,
                                                 mu_by_group,
                                                 sigma_by_group,
                                                 shared,
                                                 feature) {
  formula <- .scdesign3_compact_formula(formula)
  sigma_formula <- .scdesign3_compact_formula(sigma_formula)
  names(mu_by_group) <- shared$group_signatures
  names(sigma_by_group) <- shared$group_signatures
  fitted <- mu_by_group[shared$group_signatures[shared$train_group]]
  sigma_fitted <- sigma_by_group[shared$group_signatures[shared$train_group]]
  names(fitted) <- names(y)
  names(sigma_fitted) <- names(y)
  ll <- sum(stats::dnbinom(
    y,
    size = 1 / pmax(sigma_fitted, .Machine$double.eps),
    mu = pmax(fitted, .Machine$double.eps),
    log = TRUE
  ))
  df <- 2L * length(mu_by_group)

  structure(
    list(
      formula = formula,
      terms = .scdesign3_compact_terms(stats::terms(formula, data = shared$data)),
      mu.formula = formula,
      sigma.formula = sigma_formula,
      family = .scdesign3_categorical_nbi_family(),
      family_key = "nb",
      y = y,
      mu_by_group = mu_by_group,
      sigma_by_group = sigma_by_group,
      logLik = ll,
      df = df,
      df.residual = max(length(y) - df, 0L),
      nobs = length(y),
      shared = shared,
      backend = "scDesign3::categorical_nbi_group_mle",
      feature = feature
    ),
    class = "scdesign3_categorical_nbi"
  )
}

.scdesign3_categorical_nbi_family <- function() {
  list(
    family = "NBI",
    link = "log",
    linkfun = log,
    linkinv = function(eta) pmax(exp(eta), .Machine$double.eps)
  )
}

.scdesign3_is_categorical_nbi <- function(fit) {
  inherits(fit, "scdesign3_categorical_nbi")
}

.scdesign3_is_categorical_dist <- function(fit) {
  inherits(fit, "scdesign3_categorical_dist") || .scdesign3_is_categorical_nbi(fit)
}

.scdesign3_predict_categorical_dist <- function(object,
                                                newdata = NULL,
                                                what = c("mu", "sigma", "nu", "size")) {
  what <- match.arg(what)
  if (.scdesign3_is_categorical_nbi(object) && !inherits(object, "scdesign3_categorical_dist")) {
    return(.scdesign3_predict_categorical_nbi(object, newdata = newdata, what = what))
  }

  if (is.null(newdata)) {
    idx <- object$shared$train_group
    values <- switch(
      what,
      mu = object$mu_by_group,
      sigma = object$sigma_by_group,
      nu = object$nu_by_group,
      size = 1 / pmax(object$sigma_by_group, .Machine$double.eps)
    )
    if (is.null(values)) {
      out <- rep(if (identical(what, "nu")) 0 else NA_real_, length(idx))
    } else {
      out <- as.numeric(values[idx])
    }
    names(out) <- names(object$y)
    return(out)
  }

  x <- .scdesign3_categorical_predict_matrix(object, newdata)
  design_signatures <- .scdesign3_model_matrix_signature(x)
  idx <- match(design_signatures, object$shared$group_signatures)
  if (anyNA(idx)) {
    stop("New data contains categorical design rows not observed during fitting.", call. = FALSE)
  }
  values <- switch(
    what,
    mu = object$mu_by_group,
    sigma = object$sigma_by_group,
    nu = object$nu_by_group,
    size = 1 / pmax(object$sigma_by_group, .Machine$double.eps)
  )
  if (is.null(values)) {
    out <- rep(if (identical(what, "nu")) 0 else NA_real_, length(idx))
  } else {
    out <- as.numeric(values[idx])
  }
  names(out) <- rownames(newdata)
  out
}

.scdesign3_predict_categorical_nbi <- function(object,
                                               newdata = NULL,
                                               what = c("mu", "sigma", "size")) {
  what <- match.arg(what)
  if (is.null(newdata)) {
    idx <- object$shared$train_group
    values <- switch(
      what,
      mu = object$mu_by_group,
      sigma = object$sigma_by_group,
      size = 1 / pmax(object$sigma_by_group, .Machine$double.eps)
    )
    out <- as.numeric(values[idx])
    names(out) <- names(object$y)
    return(out)
  }

  x <- .scdesign3_categorical_predict_matrix(object, newdata)
  design_signatures <- .scdesign3_model_matrix_signature(x)
  idx <- match(design_signatures, object$shared$group_signatures)
  if (anyNA(idx)) {
    stop("New data contains categorical design rows not observed during fitting.", call. = FALSE)
  }
  values <- switch(
    what,
    mu = object$mu_by_group,
    sigma = object$sigma_by_group,
    size = 1 / pmax(object$sigma_by_group, .Machine$double.eps)
  )
  out <- as.numeric(values[idx])
  names(out) <- rownames(newdata)
  out
}

.scdesign3_categorical_predict_matrix <- function(object, newdata) {
  newdata <- as.data.frame(newdata)
  mf <- stats::model.frame(
    object$shared$rhs_terms,
    data = newdata,
    na.action = stats::na.pass,
    xlev = object$shared$xlevels
  )
  x <- stats::model.matrix(
    object$shared$rhs_terms,
    data = mf,
    contrasts.arg = object$shared$contrasts
  )
  .scdesign3_align_matrix_columns(x, object$shared$design_names)
}

predict.scdesign3_categorical_dist <- function(object,
                                               newdata = NULL,
                                               type = c("response", "link"),
                                               what = c("mu", "sigma", "nu"),
                                               ...) {
  type <- match.arg(type)
  what <- match.arg(what)
  out <- .scdesign3_predict_categorical_dist(object, newdata = newdata, what = what)
  if (identical(type, "link")) {
    if (identical(what, "mu") && identical(object$family_key, "gaussian")) {
      return(out)
    }
    if (identical(what, "mu") && identical(object$family_key, "binomial")) {
      return(stats::qlogis(pmin(pmax(out, .Machine$double.eps), 1 - .Machine$double.eps)))
    }
    if (identical(what, "nu") ||
        (identical(what, "sigma") && identical(object$family_key, "zip"))) {
      return(stats::qlogis(pmin(pmax(out, .Machine$double.eps), 1 - .Machine$double.eps)))
    }
    out <- log(pmax(out, .Machine$double.eps))
  }
  out
}

predict.scdesign3_categorical_nbi <- function(object,
                                              newdata = NULL,
                                              type = c("response", "link"),
                                              what = c("mu", "sigma"),
                                              ...) {
  type <- match.arg(type)
  what <- match.arg(what)
  out <- .scdesign3_predict_categorical_nbi(object, newdata = newdata, what = what)
  if (identical(type, "link")) {
    out <- log(pmax(out, .Machine$double.eps))
  }
  out
}

family.scdesign3_categorical_dist <- function(object, ...) {
  object$family
}

logLik.scdesign3_categorical_dist <- function(object, ...) {
  out <- object$logLik
  attr(out, "df") <- object$df
  attr(out, "nobs") <- object$nobs
  class(out) <- "logLik"
  out
}

nobs.scdesign3_categorical_dist <- function(object, ...) {
  object$nobs
}

family.scdesign3_categorical_nbi <- function(object, ...) {
  object$family
}

logLik.scdesign3_categorical_nbi <- function(object, ...) {
  out <- object$logLik
  attr(out, "df") <- object$df
  attr(out, "nobs") <- object$nobs
  class(out) <- "logLik"
  out
}

nobs.scdesign3_categorical_nbi <- function(object, ...) {
  object$nobs
}

.scdesign3_scglm_family <- function(family_key, theta = NULL) {
  switch(
    family_key,
    binomial = stats::binomial(),
    poisson = stats::poisson(),
    gaussian = stats::gaussian(),
    nb = {
      if (is.null(theta) || !is.finite(theta)) {
        mgcv::nb()
      } else {
        mgcv::nb(theta = theta)
      }
    },
    stop("Unsupported family for scGLM backend: ", family_key, call. = FALSE)
  )
}

.scdesign3_scglm_method <- function(method, family_key, formula, data) {
  method <- match.arg(
    method,
    c("auto", "categorical_closed_form", "categorical_irls", "irls", "newton_stein")
  )
  if (!identical(method, "auto")) {
    return(method)
  }
  if (identical(family_key, "nb")) {
    if (.scdesign3_formula_is_intercept_only_rhs(formula)) {
      return("newton_stein")
    }
    if (.scdesign3_formula_is_all_categorical(formula, data)) {
      return("categorical_closed_form")
    }
    return("unsupported")
  }
  if (.scdesign3_formula_is_all_categorical(formula, data)) {
    return("categorical_closed_form")
  }
  "irls"
}

.scdesign3_formula_is_intercept_only_rhs <- function(formula) {
  rhs_terms <- stats::terms(formula)
  length(attr(rhs_terms, "term.labels")) == 0L && attr(rhs_terms, "intercept") == 1L
}

.scdesign3_formula_is_all_categorical <- function(formula, data) {
  if (grepl("s\\(|te\\(", paste(deparse(formula), collapse = " "))) {
    return(FALSE)
  }
  rhs_terms <- stats::terms(formula, data = data)
  rhs_frame <- stats::model.frame(rhs_terms, data = data, na.action = stats::na.omit)
  ncol(rhs_frame) > 0L && all(vapply(rhs_frame, function(x) {
    is.factor(x) || is.character(x) || is.logical(x) ||
      (is.numeric(x) && all(x == floor(x), na.rm = TRUE))
  }, logical(1)))
}

.scdesign3_formula_has_smooth <- function(formula) {
  text <- paste(deparse(stats::as.formula(paste0("~", as.character(formula)[1L]))), collapse = " ")
  grepl("\\b(s|te|ti|t2)\\s*\\(", text, perl = TRUE)
}

.scdesign3_is_intercept_only_formula <- function(x) {
  x <- gsub("\\s+", "", as.character(x)[1])
  x %in% c("1", "~1")
}

.scdesign3_scglm_training_x <- function(fit) {
  if (!is.null(fit$x)) {
    return(fit$x)
  }
  if (!is.null(fit$x_obs)) {
    return(fit$x_obs)
  }
  stop("Could not extract the scGLM training design matrix.", call. = FALSE)
}

.scdesign3_scglm_basis_object <- function(formula, data) {
  if (!grepl("s\\(|te\\(", paste(deparse(formula), collapse = " "))) {
    return(NULL)
  }
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    return(NULL)
  }
  basis_data <- as.data.frame(data)
  # A constant Gaussian response can make mgcv's REML fit singular for some
  # smooth bases. The response is only used to create a reusable lpmatrix.
  basis_data$gene <- seq_len(nrow(basis_data)) / max(nrow(basis_data), 1L)
  tryCatch(
    mgcv::gam(formula = formula, data = basis_data, family = stats::gaussian(), method = "REML"),
    error = function(e) NULL
  )
}

.scdesign3_new_scglm_model <- function(formula,
                                       rhs_formula,
                                       family_key,
                                       coefficients,
                                       x,
                                       y,
                                       fitted,
                                       eta,
                                       theta,
                                       data,
                                       shared,
                                       backend,
                                       feature) {
  coef <- as.numeric(coefficients)
  names(coef) <- rownames(coefficients)
  if (is.null(names(coef))) {
    names(coef) <- colnames(x)
  }
  sig2 <- if (identical(family_key, "gaussian")) {
    sum((y - fitted)^2) / max(length(y) - length(coef), 1L)
  } else {
    1
  }
  ll <- .scdesign3_scglm_loglik(y = y, mu = fitted, family_key = family_key, theta = theta, sig2 = sig2)
  df <- length(coef) + as.integer(identical(family_key, "nb")) + as.integer(identical(family_key, "gaussian"))
  fam <- .scdesign3_scglm_model_family(family_key, theta = theta)
  terms_obj <- stats::terms(formula, data = data)

  structure(
    list(
      formula = formula,
      terms = terms_obj,
      rhs_formula = rhs_formula,
      family = fam,
      family_key = family_key,
      coefficients = coef,
      fitted.values = fitted,
      linear.predictors = eta,
      y = y,
      residuals = y - fitted,
      sig2 = sig2,
      theta = theta,
      logLik = ll,
      df = df,
      edf = rep(1, length(coef)),
      edf1 = rep(1, length(coef)),
      edf2 = rep(1, length(coef)),
      df.residual = max(length(y) - df, 0L),
      nobs = length(y),
      shared = shared,
      backend = backend,
      feature = feature
    ),
    class = "scdesign3_scglm"
  )
}

.scdesign3_scglm_loglik <- function(y, mu, family_key, theta = NA_real_, sig2 = 1) {
  mu <- pmax(as.numeric(mu), .Machine$double.eps)
  switch(
    family_key,
    binomial = sum(stats::dbinom(y, size = 1, prob = pmin(pmax(mu, 1e-12), 1 - 1e-12), log = TRUE)),
    poisson = sum(stats::dpois(y, lambda = mu, log = TRUE)),
    gaussian = sum(stats::dnorm(y, mean = mu, sd = sqrt(max(sig2, .Machine$double.eps)), log = TRUE)),
    nb = sum(stats::dnbinom(y, size = max(theta, .Machine$double.eps), mu = mu, log = TRUE)),
    stop("Unsupported family for scGLM log-likelihood.", call. = FALSE)
  )
}

.scdesign3_scglm_model_family <- function(family_key, theta = NA_real_) {
  switch(
    family_key,
    binomial = stats::binomial(),
    poisson = stats::poisson(),
    gaussian = stats::gaussian(),
    nb = {
      th <- as.numeric(theta)[1L]
      list(
        family = paste0("Negative Binomial(", signif(th, 6), ")"),
        link = "log",
        linkinv = function(eta) pmax(exp(eta), .Machine$double.eps),
        getTheta = function(trans = FALSE) th
      )
    },
    stop("Unsupported family for scGLM model object: ", family_key, call. = FALSE)
  )
}

.scdesign3_filtered_marginal <- function(gene, trace) {
  warning_log <- list(list(
    function_name = "fit_marginal",
    type = "warning",
    message = paste0(gene, "is expressed in too few cells.")
  ))
  if (isTRUE(trace)) {
    return(list(fit = NA, warning = warning_log, time = c(NA_real_, NA_real_), removed_cell = NA))
  }
  list(fit = NA, removed_cell = NA)
}

predict.scdesign3_scglm <- function(object,
                                    newdata = NULL,
                                    type = c("link", "response"),
                                    ...) {
  type <- match.arg(type)
  eta <- if (is.null(newdata)) {
    object$linear.predictors
  } else {
    x <- .scdesign3_scglm_predict_matrix(object, newdata)
    as.vector(x %*% object$coefficients)
  }
  if (identical(type, "response")) {
    out <- object$family$linkinv(eta)
  } else {
    out <- eta
  }
  if (is.null(names(out))) {
    names(out) <- if (is.null(newdata)) names(object$fitted.values) else rownames(newdata)
  }
  out
}

.scdesign3_scglm_predict_matrix <- function(object, newdata) {
  newdata <- as.data.frame(newdata)
  shared <- .scdesign3_scglm_shared(object)
  if (!is.null(shared$basis_object)) {
    basis_data <- newdata
    basis_data$gene <- 0
    lp <- stats::predict(shared$basis_object, newdata = basis_data, type = "lpmatrix")
    return(.scdesign3_align_matrix_columns(lp, shared$design_names))
  }
  mf <- stats::model.frame(
    shared$rhs_terms,
    data = newdata,
    na.action = stats::na.pass,
    xlev = shared$xlevels
  )
  x <- stats::model.matrix(shared$rhs_terms, data = mf, contrasts.arg = shared$contrasts)
  .scdesign3_align_matrix_columns(x, shared$design_names)
}

.scdesign3_align_matrix_columns <- function(x, target_names) {
  out <- matrix(0, nrow = nrow(x), ncol = length(target_names))
  colnames(out) <- target_names
  rownames(out) <- rownames(x)
  idx <- match(target_names, colnames(x))
  keep <- !is.na(idx)
  out[, keep] <- x[, idx[keep], drop = FALSE]
  out
}

family.scdesign3_scglm <- function(object, ...) {
  object$family
}

logLik.scdesign3_scglm <- function(object, ...) {
  out <- object$logLik
  attr(out, "df") <- object$df
  attr(out, "nobs") <- object$nobs
  class(out) <- "logLik"
  out
}

nobs.scdesign3_scglm <- function(object, ...) {
  object$nobs
}

.scdesign3_scglm_shared <- function(object) {
  if (!is.null(object$shared)) {
    return(object$shared)
  }
  object
}

.scdesign3_xlevels <- function(terms, data) {
  vars <- all.vars(terms)
  vars <- intersect(vars, names(data))
  out <- lapply(data[vars], function(x) {
    if (is.factor(x)) {
      levels(x)
    } else {
      NULL
    }
  })
  out[!vapply(out, is.null, logical(1))]
}
