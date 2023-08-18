find_genes <- function(myRCTD, cell_types, MEAN_VAL_THRESH=5, TOT_THRESH=100) {
  mean_val <- myRCTD@de_results$gene_fits$mean_val # Cell type mean gene expression from CSIDE
  mean_val <- mean_val[,cell_types]
  # find genes with well-behaved estimates i.e. low variance
  tot <- rowSums(myRCTD@spatialRNA@maternalCounts) + rowSums(myRCTD@spatialRNA@paternalCounts)
  tot <- tot[rownames(mean_val)]
  genes <- rownames(mean_val)[which(rowVars(mean_val)<MEAN_VAL_THRESH & tot>TOT_THRESH)]
  return(genes)
}


simulate_data <- function(myRCTD, gene, cell_type_with_pattern,
                          cell_types, df=5, seed=1, add_count = 20,
                          phi_j=0.3, doublet_mode=F, npix = 1e4,
                          my_beta = NULL, X2=NULL) {
  set.seed(seed)
  puck <- myRCTD@originalSpatialRNA
  if (is.null(X2)) {
    X2 <- build.designmatrix.nonparam(myRCTD, df = df)
  }
  barcodes <- rownames(X2)
  if (doublet_mode) {
    nonreject <- barcodes[which(myRCTD@results$results_df$spot_class != 'reject' & myRCTD@results$results_df$first_type %in% cell_types)]
    barcodes <- sample(nonreject, npix, replace=F)
    X2 <- X2[barcodes,]
  }
  coords <- myRCTD@spatialRNA@coords[barcodes,]
  coef <- rnorm(df, 0, 2)
  p <- expit(X2%*%coef)
  if (doublet_mode) {
    if (is.null(my_beta)) {
      my_beta <- get_beta_doublet(barcodes, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
    }
    thresh <- 0.999
  } else {
    weights <- myRCTD@results$weights
    weights <- weights[,cell_types]
    my_beta <- as.matrix(sweep(weights, 1, rowSums(weights), '/')) # Cell type weights from full mode RCTD
    thresh <- 0.8
  }
  mean_val <- myRCTD@de_results$gene_fits$mean_val # Cell type mean gene expression from CSIDE
  mean_val <- matrix(pmax(-10,pmin(10,mean_val)),
                     nrow = nrow(mean_val),
                     ncol=ncol(mean_val),
                     dimnames = dimnames(mean_val))# numerical stability
  mean_val <- mean_val[,cell_types]
  alphas <- sweep(my_beta[barcodes,cell_types], 2, exp(mean_val[gene,]), '*')
  alphas <- sweep(alphas, 1, rowSums(alphas), '/')
  if (doublet_mode) {
    nanidx <- which(is.na(rowSums(alphas)))
    if (length(nanidx) > 0) {
      barcodes <- barcodes[-nanidx]
    }
    X2 <- X2[barcodes,]
    my_beta <- my_beta[barcodes,]
    alphas <- alphas[barcodes,]
  }
  ## Compute p_ij
  p_ij <- numeric(nrow(alphas))
  for (k in 1:ncol(alphas)) {
    if (colnames(alphas)[k] == cell_type_with_pattern) {
      p_ij <- p_ij + alphas[,k]*expit(X2%*%coef)
    } else {
      p_ij <- p_ij + alphas[,k]*expit(0) # expit(0) = 0.5
    }
  }
  ## Convert to Beta(a,b)
  a <- p_ij*(1-phi_j)/phi_j
  b <- (1-phi_j)*(1-p_ij)/phi_j
  lambd_ij <- rbeta(length(a), a, b)
  ## Binomial sampling
  N_ij <- myRCTD@spatialRNA@maternalCounts[gene,barcodes]+myRCTD@spatialRNA@paternalCounts[gene,barcodes]+add_count
  y_ij <- rbinom(length(a), size = N_ij, prob = lambd_ij)
  return(list(design_mat=X2, coef=coef,y=y_ij, N=N_ij, barcodes=barcodes, my_beta=my_beta))
}

estimate_single_gene <- function(myRCTD, gene, barcodes, y, N, my_beta=NULL, X2=NULL, df=5, doublet_mode = F) {
  if (doublet_mode) {
    myRCTD@originalSpatialRNA@maternalCounts[gene,] <- 0
    myRCTD@originalSpatialRNA@paternalCounts[gene,] <- 0
  }
  myRCTD@originalSpatialRNA@maternalCounts[gene,barcodes] <- y
  myRCTD@originalSpatialRNA@paternalCounts[gene,barcodes] <- N-y
  fit_genes_subset = gene
  barcodes = NULL;
  cell_type_threshold = 128; gene_threshold = 5e-5;
  weight_threshold = NULL; sigma_gene = T;
  PRECISION.THRESHOLD = 0.05; cell_types_present = NULL; fdr = .01; test_genes_sig = T;
  logs=F; test_error = F; spase = T; remove_genes = NULL

  if (!is.null(X2)) {
    X2 <- X2
  } else {
    X2 <- build.designmatrix.nonparam(myRCTD, barcodes = barcodes, df = df)
  }
  region_thresh <- cell_type_threshold / 4
  barcodes <- rownames(X2)
  coords <- myRCTD@spatialRNA@coords[barcodes,]
  medx <- median(coords$x); medy <- median(coords$y)
  r1 <- barcodes[coords$x < medx & coords$y < medy]
  cell_type_filter <- aggregate_cell_types(myRCTD, r1, doublet_mode = doublet_mode) >= region_thresh
  r2 <- barcodes[coords$x < medx & coords$y > medy]
  cell_type_filter <- cell_type_filter & (aggregate_cell_types(myRCTD, r2, doublet_mode = doublet_mode) >= region_thresh)
  r3 <- barcodes[coords$x > medx & coords$y > medy]
  cell_type_filter <- cell_type_filter & (aggregate_cell_types(myRCTD, r3, doublet_mode = doublet_mode) >= region_thresh)
  r4 <- barcodes[coords$x > medx & coords$y > medy]
  cell_type_filter <- cell_type_filter & (aggregate_cell_types(myRCTD, r4, doublet_mode = doublet_mode) >= region_thresh)
  cell_type_count <- aggregate_cell_types(myRCTD, barcodes, doublet_mode = doublet_mode)

  test_mode = 'individual'
  params_to_test = 2:df; normalize_expr = F;
  cell_type_specific <- NULL

  X <- check_designmatrix(X2, 'run.CSIDE', require_2d = TRUE)
  if(is.null(cell_type_specific))
    cell_type_specific <- !logical(dim(X)[2])
  check_cell_type_specific(cell_type_specific, dim(X)[2])
  X1 <- X[,!cell_type_specific];
  if(any(!cell_type_specific)) {
    X2 <- X[,cell_type_specific]
  } else {
    X2 <- X
  }
  log_fc_thresh = 0; test_error = FALSE; fdr_method = 'BH';


  if(gene_threshold == .01 || fdr == 0.25 || cell_type_threshold == 10 ||
     (!is.null(weight_threshold) && weight_threshold == 0.1))
    warning('run.CSIDE.general: some parameters are set to the CSIDE vignette values, which are intended for testing but not proper execution. For more accurate results, consider using the default parameters to this function.')
  if(doublet_mode && myRCTD@config$RCTDmode != 'doublet')
    stop('run.CSIDE.general: attempted to run CSIDE in doublet mode, but RCTD was not run in doublet mode. Please run CSIDE in full mode (doublet_mode = F) or run RCTD in doublet mode.')
  if(!any("cell_types_assigned" %in% names(myRCTD@internal_vars)) || !myRCTD@internal_vars$cell_types_assigned)
    stop('run.CSIDE.general: cannot run CSIDE unless cell types have been assigned. If cell types have been assigned, you may run "myRCTD <- set_cell_types_assigned(myRCTD)".')
  if((myRCTD@config$RCTDmode != 'multi') && (length(setdiff(barcodes,rownames(myRCTD@results$weights))) > 0)) {
    warning('run.CSIDE.general: some elements of barcodes do not appear in myRCTD object (myRCTD@results$weights), but they are required to be a subset. Downsampling barcodes to the intersection of the two sets.')
    barcodes <- intersect(barcodes,rownames(myRCTD@results$weights))
  }
  cell_type_info <- myRCTD@cell_type_info$info
  if(doublet_mode) {
    if (is.null(my_beta)) {
      my_beta <- get_beta_doublet(barcodes, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
    }
    thresh <- 0.999
  } else if(myRCTD@config$RCTDmode == "multi") {
    my_beta <- get_beta_multi(barcodes, cell_type_info[[2]], myRCTD@results, myRCTD@spatialRNA@coords)
    thresh <- 0.999
  } else {
    my_beta <- as.matrix(sweep(myRCTD@results$weights, 1, rowSums(myRCTD@results$weights), '/'))
    thresh <- 0.8
  }
  if(!is.null(weight_threshold))
    thresh <- weight_threshold
  thresh <- 0.5
  cell_type_threshold <- 2
  cell_types <- choose_cell_types(myRCTD, barcodes, doublet_mode, cell_type_threshold, cell_types,
                                  my_beta, thresh, cell_type_filter)
  if(length(cell_types) < 1)
    stop('run.CSIDE.general: zero cell types remain. Cannot run CSIDE with zero cell types.')
  X1 <- check_designmatrix(X1, 'run.CSIDE.general')
  X2 <- check_designmatrix(X2, 'run.CSIDE.general', require_2d = TRUE)
  if(!(test_mode %in% c('individual', 'categorical')))
    stop(c('run.CSIDE.general: not valid test_mode = ',test_mode,'. Please set test_mode = "categorical" or "individual".'))
  if(is.null(params_to_test)) {
    if(test_mode == 'individual') {
      params_to_test <- min(2, dim(X2)[2])
    } else {
      params_to_test <- 1:dim(X2)[2]
    }
  }

  if(normalize_expr && (test_mode != 'individual' || length(params_to_test) > 1))
    stop('run.CSIDE.general: Setting normalize_expr = TRUE is only valid for testing single parameters with test_mode = individual')
  if(any(!(params_to_test %in% 1:dim(X2)[2])))
    stop(c('run.CSIDE.general: params_to_test must be a vector of integers from 1 to dim(X2)[2] = ', dim(X2)[2],
           'please make sure that tested parameters are in the required range.'))
  if(test_mode == 'categorical' && any(!(X2[,params_to_test] %in% c(0,1))))
    stop(c('run.CSIDE.general: for test_mode = categorical, colums params_to_test, ',params_to_test,', must have values 0 or 1.'))
  if(is.null(cell_types_present))
    cell_types_present <- cell_types
  if(any(!(barcodes %in% rownames(X1))) || any(!(barcodes %in% rownames(X2))))
    stop('run.CSIDE.general: some barcodes do not appear in the rownames of X1 or X2.')
  puck = myRCTD@originalSpatialRNA
  gene_list_tot <- gene
  if(length(gene_list_tot) == 0)
    stop('run.CSIDE.general: no genes past threshold. Please consider lowering gene_threshold.')
  if (spase) {
    nUMI <- colSums(puck@maternalCounts) + colSums(puck@paternalCounts)
    names(nUMI) <- colnames(puck@maternalCounts)
    puck@nUMI <- nUMI
    nUMI <- nUMI[barcodes]
  } else {
    nUMI <- puck@nUMI[barcodes]
  }

  res <- filter_barcodes_cell_types(barcodes, cell_types, my_beta, thresh = thresh)
  if(test_error)
    return(myRCTD)
  barcodes <- res$barcodes; my_beta <- res$my_beta
  if(!spase) {
    sigma_init <- as.character(100*myRCTD@internal_vars$sigma)
    if(sigma_gene) {
      set_global_Q_all()
      sigma_set <- sigma_init
      set_likelihood_vars(Q_mat_all[[sigma_init]], X_vals, sigma = sigma_set)
    } else {
      set_likelihood_vars_sigma(sigma_init)
    }
  } else {
    sigma_init <- '1'
  }

  X1 <- X1[barcodes, , drop = FALSE]
  X2 <- X2[barcodes, , drop = FALSE]
  puck <- restrict_puck(puck, barcodes)

  mean_val <- myRCTD@de_results$gene_fits$mean_val[,cell_types]
  mean_val <- matrix(pmax(-10,pmin(10,mean_val)),
                     nrow = nrow(mean_val),
                     ncol=ncol(mean_val),
                     dimnames = dimnames(mean_val))# numerical stability
  Y <- puck@maternalCounts[gene, barcodes]
  nUMI <- Y + puck@paternalCounts[gene,barcodes]
  my_beta_updated <- sweep(my_beta, 2, exp(mean_val[gene,]), '*')
  my_beta_updated <- sweep(my_beta_updated, 1, rowSums(my_beta_updated), '/')

  alpha1_init <- NULL; alpha2_init <- NULL
  test_mode <- 'individual'
  n.iter = 200
  MIN_CHANGE = 1e-3
  sigma_init = 1
  r = sigma_init
  MIN_ITERATIONS = 15
  PRECISION.THRESHOLD = 0.05
  res <- solveIRWLS.effects_trust(Y,X1,X2, my_beta_updated, test_mode, verbose = verbose,
                                  n.iter = n.iter, MIN_CHANGE = MIN_CHANGE,
                                  PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                                  alpha1_init = alpha1_init, alpha2_init = alpha2_init,
                                  MIN_ITERATIONS = MIN_ITERATIONS, spase = T, r = r, nUMI = nUMI)
  return(list(est = res$alpha2, var = diag(res$I)))
}

compute_metrics <- function(res, data, est_col) {
  # 1. Total sample size
  N <- sum(data$N)
  # 2. Average sample size per pixel
  aN <- mean(data$N)
  # 3. Median sample size per pixel
  mN <- median(data$N)
  # 4. MSE of coefficient point estimation
  est_coef <- res$est[,est_col]
  mse <- mean((est_coef-data$coef)^2)
  # 5. Average variance of coefficient point estimate
  idx <- (1+(est_col-1)*5):(est_col*5)
  avg_var <- mean(res$var[idx])
  # 6. Coverage (whether true coef is contained in 95% CI)
  ci_lower <- est_coef - 1.96*res$var[idx]
  ci_upper <- est_coef + 1.96*res$var[idx]
  coverage <- sum(data$coef > ci_lower & data$coef < ci_upper)
  # 7. Correlation of estimated logit(p)'s with ground truth logit(p)'s
  cor_p <- cor(data$design_mat %*% est_coef, data$design_mat %*% data$coef)
  return(c(N, aN, mN, mse, avg_var, coverage, cor_p))
}
