library(spacexr)
library(spASE)

# fixing cere 4 task 3337




# fixing cere 3 task 1726

cere3 <- readRDS('cside_spase_combined_cere_3_df_5.rds')
cere3@config$max_cores <- 1
myRCTD <- cere3
cell_types <- c('Fibroblast', 'Granule', 'MLI2','Oligodendrocytes')
df = 5; barcodes = NULL;
cell_type_threshold = 128; gene_threshold = 5e-5; doublet_mode = F;
weight_threshold = NULL; sigma_gene = T;
PRECISION.THRESHOLD = 0.05; cell_types_present = NULL; fdr = .01; test_genes_sig = T;
logs=F; test_error = F; spase = F; remove_genes = NULL

X2 <- build.designmatrix.nonparam(myRCTD, barcodes = barcodes, df = df)
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
log_fc_thresh = 0.4; test_error = FALSE; fdr_method = 'BH';


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
  my_beta <- get_beta_doublet(barcodes, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
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
cell_types <- choose_cell_types(myRCTD, barcodes, doublet_mode, cell_type_threshold, cell_types,
                                my_beta, thresh, cell_type_filter)
if(length(cell_types) < 1)
  stop('run.CSIDE.general: zero cell types remain. Cannot run CSIDE with zero cell types.')
message(paste0("run.CSIDE.general: running CSIDE with cell types ",paste(cell_types, collapse = ', ')))
X1 <- check_designmatrix(X1, 'run.CSIDE.general')
X2 <- check_designmatrix(X2, 'run.CSIDE.general', require_2d = TRUE)
if(!(test_mode %in% c('individual', 'categorical')))
  stop(c('run.CSIDE.general: not valid test_mode = ',test_mode,'. Please set test_mode = "categorical" or "individual".'))
if(is.null(params_to_test))
  if(test_mode == 'individual')
    params_to_test <- min(2, dim(X2)[2])
else
  params_to_test <- 1:dim(X2)[2]
if(normalize_expr && (test_mode != 'individual' || length(params_to_test) > 1))
  stop('run.CSIDE.general: Setting normalize_expr = TRUE is only valid for testing single parameters with test_mode = individual')
message(paste0("run.CSIDE.general: configure params_to_test = ",
               paste(paste0(params_to_test, ', ', collapse = ""))))
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
gene_list_tot <- filter_genes(puck, threshold = gene_threshold, remove_genes = remove_genes)
if(length(gene_list_tot) == 0)
  stop('run.CSIDE.general: no genes past threshold. Please consider lowering gene_threshold.')
if(length(intersect(gene_list_tot,rownames(myRCTD@cell_type_info$info[[1]]))) == 0)
  stop('run.CSIDE.general: no genes that past threshold were contained in the single cell reference. Please lower gene threshold or ensure that there is agreement between the single cell reference genes and the SpatialRNA genes.')
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

mean_val <- myRCTD@de_results$gene_fits$mean_val
mean_val <- matrix(pmax(-50,pmin(50,mean_val)),
                   nrow = nrow(mean_val),
                   ncol=ncol(mean_val),
                   dimnames = dimnames(mean_val))# numerical stability
gene <- 'Msmo1'
Y <- puck@maternalCounts[gene, barcodes]
nUMI <- Y + puck@paternalCounts[gene,barcodes]
my_beta_updated <- sweep(my_beta, 2, exp(mean_val[gene,]), '*')
my_beta_updated <- sweep(my_beta_updated, 1, rowSums(my_beta_updated), '/')
#res <- estimate_gene_wrapper(Y,X1,X2,my_beta_updated, nUMI, sigma_init, 'individual',
#                             verbose = F, n.iter = 200, MIN_CHANGE = 1e-3,
#                             sigma_gene = T, PRECISION.THRESHOLD = .05, spase = T)


alpha1_init <- NULL; alpha2_init <- NULL
test_mode <- 'individual'
n.iter = 200
MIN_CHANGE = 1e-3
sigma_init = 1
r = sigma_init
MIN_ITERATIONS = 15
PRECISION.THRESHOLD = 0.05
solveIRWLS.effects_trust(Y,X1,X2, my_beta_updated, test_mode, verbose = verbose,
                         n.iter = n.iter, MIN_CHANGE = MIN_CHANGE,
                         PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                         alpha1_init = alpha1_init, alpha2_init = alpha2_init,
                         MIN_ITERATIONS = MIN_ITERATIONS, spase = T, r = r, nUMI = nUMI)



load('cside_spase_combined_hippo_1_df_5.RData')

y = myhippo_intercept@spatialRNA@maternalCounts['Meg3',]
tot = myhippo_intercept@spatialRNA@maternalCounts['Meg3',] + myhippo_intercept@spatialRNA@paternalCounts['Meg3',]

mean(y)
mean(tot)


myhippo_intercept@spase_results$gene_fits$mean_val['Meg3',]


# starting from RCTD object


myhippo <- readRDS('rctd_combined_hippo_1.rds')
myhippo@config$max_cores <- 1
cell_types <- c('Astrocyte', 'CA1', 'CA3', 'Dentate', 'Interneuron', 'Oligodendrocyte')
myhippo_intercept <- run.CSIDE.intercept(myhippo, cell_types = cell_types, cell_type_threshold = 50, spase=F, logs=T)
myhippo_intercept <- run.CSIDE.intercept(myhippo_intercept, cell_types = cell_types, cell_type_threshold = 50, spase=T, logs=T)
save(myhippo, myhippo_intercept, file='cside_spase_combined_hippo_1_df_5.RData')



# Ptgds with intercept
doublet_mode <- T
myRCTD <- myhippo_intercept
cell_type_info <- myRCTD@cell_type_info$info
barcodes <- rownames(myRCTD@results$results_df)
if(doublet_mode) {
  my_beta <- get_beta_doublet(barcodes, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
  thresh <- 0.999
} else if(myRCTD@config$RCTDmode == "multi") {
  my_beta <- get_beta_multi(barcodes, cell_type_info[[2]], myRCTD@results, myRCTD@spatialRNA@coords)
  thresh <- 0.999
} else {
  my_beta <- as.matrix(sweep(myRCTD@results$weights, 1, rowSums(myRCTD@results$weights), '/'))
  thresh <- 0.8
}
weight_threshold <- NULL
if(!is.null(weight_threshold))
  thresh <- weight_threshold
cell_type_threshold <- 50
cell_type_filter <- NULL
cell_types <- myRCTD@internal_vars_de$cell_types
cell_types <- choose_cell_types(myRCTD, barcodes, doublet_mode, cell_type_threshold, cell_types,
                                my_beta, thresh, cell_type_filter)
if(length(cell_types) < 1)
  stop('run.CSIDE.general: zero cell types remain. Cannot run CSIDE with zero cell types.')
message(paste0("run.CSIDE.general: running CSIDE with cell types ",paste(cell_types, collapse = ', ')))
X2 <- build.designmatrix.intercept(myRCTD, barcodes = barcodes)
X1 <- X2[,0:0]
X1 <- check_designmatrix(X1, 'run.CSIDE.general')
X2 <- check_designmatrix(X2, 'run.CSIDE.general', require_2d = TRUE)
cell_types_present <- cell_types
puck = myRCTD@originalSpatialRNA
nUMI <- colSums(puck@maternalCounts) + colSums(puck@paternalCounts)
names(nUMI) <- colnames(puck@maternalCounts)
nUMI <- nUMI[barcodes]
res <- filter_barcodes_cell_types(barcodes, cell_types, my_beta, thresh = thresh)
barcodes <- res$barcodes; my_beta <- res$my_beta
sigma_init <- '1'

X1 <- X1[barcodes, , drop = FALSE]
X2 <- X2[barcodes, , drop = FALSE]
#,my_beta, nUMI[barcodes], gene_list_tot,
#                                cell_types, restrict_puck(puck, barcodes), barcodes, sigma_init,
#                                test_mode, myRCTD, numCores = myRCTD@config$max_cores, sigma_gene = sigma_gene,
#                                PRECISION.THRESHOLD = PRECISION.THRESHOLD, params_to_test = params_to_test,
#                               logs=logs, spase = spase)


mean_val <- myRCTD@de_results$gene_fits$mean_val
gene <- 'Celf5'
Y <- puck@maternalCounts[gene, barcodes]
nUMI <- puck@maternalCounts[gene, barcodes] + puck@paternalCounts[gene,barcodes]
my_beta_updated <- sweep(my_beta, 2, exp(mean_val[gene,]), '*')
my_beta_updated <- sweep(my_beta_updated, 1, rowSums(my_beta_updated), '/')
res <- estimate_gene_wrapper(Y,X1,X2,my_beta_updated, nUMI, sigma_init, 'individual',
                             verbose = F, n.iter = 200, MIN_CHANGE = 1e-3,
                             sigma_gene = T, PRECISION.THRESHOLD = .05, spase = T)

hippo1 <- readRDS('cside_spase_combined_hippo_1_df_15.rds')
myRCTD <- hippo1
# Cpne3 with 15 df
doublet_mode <- T
cell_type_info <- myRCTD@cell_type_info$info
barcodes <- rownames(myRCTD@results$results_df)
if(doublet_mode) {
  my_beta <- get_beta_doublet(barcodes, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
  thresh <- 0.999
} else if(myRCTD@config$RCTDmode == "multi") {
  my_beta <- get_beta_multi(barcodes, cell_type_info[[2]], myRCTD@results, myRCTD@spatialRNA@coords)
  thresh <- 0.999
} else {
  my_beta <- as.matrix(sweep(myRCTD@results$weights, 1, rowSums(myRCTD@results$weights), '/'))
  thresh <- 0.8
}
weight_threshold <- NULL
if(!is.null(weight_threshold))
  thresh <- weight_threshold
cell_type_threshold <- 50
cell_type_filter <- NULL
cell_types <- myRCTD@internal_vars_de$cell_types
cell_types <- choose_cell_types(myRCTD, barcodes, doublet_mode, cell_type_threshold, cell_types,
                                my_beta, thresh, cell_type_filter)
if(length(cell_types) < 1)
  stop('run.CSIDE.general: zero cell types remain. Cannot run CSIDE with zero cell types.')
message(paste0("run.CSIDE.general: running CSIDE with cell types ",paste(cell_types, collapse = ', ')))
X2 <- build.designmatrix.nonparam(myRCTD, barcodes = barcodes, df = 15)

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

cell_types_present <- cell_types
puck = myRCTD@originalSpatialRNA
nUMI <- colSums(puck@maternalCounts) + colSums(puck@paternalCounts)
names(nUMI) <- colnames(puck@maternalCounts)
nUMI <- nUMI[barcodes]
res <- filter_barcodes_cell_types(barcodes, cell_types, my_beta, thresh = thresh)
barcodes <- res$barcodes; my_beta <- res$my_beta
sigma_init <- '1'
X1 <- X2[,0:0]
X1 <- X1[barcodes, , drop = FALSE]
X2 <- X2[barcodes, , drop = FALSE]
#,my_beta, nUMI[barcodes], gene_list_tot,
#                                cell_types, restrict_puck(puck, barcodes), barcodes, sigma_init,
#                                test_mode, myRCTD, numCores = myRCTD@config$max_cores, sigma_gene = sigma_gene,
#                                PRECISION.THRESHOLD = PRECISION.THRESHOLD, params_to_test = params_to_test,
#                               logs=logs, spase = spase)


mean_val <- myRCTD@de_results$gene_fits$mean_val
gene <- 'Cetn3'
Y <- puck@maternalCounts[gene, barcodes]
nUMI <- puck@maternalCounts[gene, barcodes] + puck@paternalCounts[gene,barcodes]
my_beta_updated <- sweep(my_beta, 2, exp(mean_val[gene,]), '*')
my_beta_updated <- sweep(my_beta_updated, 1, rowSums(my_beta_updated), '/')
#res <- estimate_gene_wrapper(Y,X1,X2,my_beta_updated, nUMI, sigma_init, 'individual',
#                             verbose = F, n.iter = 200, MIN_CHANGE = 1e-3,
#                             sigma_gene = T, PRECISION.THRESHOLD = .05, spase = T)


alpha1_init <- NULL; alpha2_init <- NULL
test_mode <- 'individual'
n.iter = 200
MIN_CHANGE = 1e-3
sigma_init = 1
r = sigma_init
MIN_ITERATIONS = 15
PRECISION.THRESHOLD = 0.05
solveIRWLS.effects_trust(Y,X1,X2, my_beta, test_mode, verbose = verbose,
                         n.iter = n.iter, MIN_CHANGE = MIN_CHANGE,
                         PRECISION.THRESHOLD = PRECISION.THRESHOLD,
                         alpha1_init = alpha1_init, alpha2_init = alpha2_init,
                         MIN_ITERATIONS = MIN_ITERATIONS, spase = T, r = r, nUMI = nUMI)
##############

devtools::load_all()

load('cside_spase_combined_hippo_1_df_5.RData')

myhippo_intercept@spase_results$gene_fits$mean_val['Meg3',]
myRCTD <- myhippo_intercept
puck <- myRCTD@originalSpatialRNA
maternal
X2 <- build.designmatrix.intercept(myhippo_intercept, barcodes = NULL)
barcodes <- rownames(X2)
gene_threshold <- 5e-5
remove_genes <- NULL

cell_type = 'Oligodendrocyte'
cell_type_info <- myRCTD@cell_type_info$info
cell_types = myhippo_intercept@internal_vars_de$cell_types
my_beta <- get_beta_doublet(barcodes, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
thresh <- 0.999
nUMI <- colSums(puck@maternalCounts) + colSums(puck@paternalCounts)
names(nUMI) <- colnames(puck@maternalCounts)
nUMI <- nUMI[barcodes]


res <- filter_barcodes_cell_types(barcodes, cell_types, my_beta, thresh = 0.999)
barcodes <- res$barcodes
my_beta <- res$my_beta

gene_fits = myhippo_intercept@spase_results$gene_fits
gene_list_tot <- filter_genes(puck, threshold = gene_threshold, remove_genes = remove_genes)
cti_renorm <- get_norm_ref(
  puck,
  myRCTD@cell_type_info$info[[1]],
  intersect(gene_list_tot,rownames(myRCTD@cell_type_info$info[[1]])),
  myRCTD@internal_vars$proportions
)

gene_list_type <- get_gene_list_type(
  my_beta = my_beta,
  barcodes = barcodes,
  cell_type = cell_type,
  nUMI = nUMI,
  gene_list_type = gene_list_tot,
  cti_renorm = cti_renorm,
  cell_types_present = cell_types,
  gene_fits = gene_fits,
  test_mode = 'individual'
)

both_genes <- find_sig_genes_individual(
  cell_type = 'Oligodendrocyte',
  cell_types = cell_types,
  gene_fits = gene_fits,
  gene_list_type = gene_list_type,
  X2 = X2,
  params_to_test = 1,
  fdr = 0.01,
  p_thresh = 1,
  log_fc_thresh = 0.4,
  normalize_expr = F,
  fdr_method = 'BH',
  se_thresh = 1
)


