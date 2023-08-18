# first argument is the chunk of params to run (1-18)
library(spacexr)
library(spASE)
devtools::load_all()

args <- commandArgs(trailingOnly=T)

source('../R/spase_sim_utils.R')


myRCTD <- readRDS('cside_spase_combined_cere_4_df_5.rds')
params <- read.delim('simulations_cere4_params.tsv')
df <- 5
cell_types <- c('Fibroblast', 'Granule', 'MLI2', 'Purkinje', 'Bergmann', 'Oligodendrocytes')

chunk <- as.integer(args[1])
idx <- (1:1000)+(chunk-1)*1000

X2 <- build.designmatrix.nonparam(myRCTD, barcodes = NULL, df = df)
barcodes <- rownames(X2)
cell_type_info <- myRCTD@cell_type_info$info
my_beta <- get_beta_doublet(barcodes, cell_type_info[[2]], myRCTD@results$results_df, myRCTD@results$weights_doublet)
for (i in idx) {
  message((paste(Sys.time(), i)))
  ct_idx <- which(cell_types == params$cell_type[i])
  data <- simulate_data(myRCTD, params$gene[i], cell_types[ct_idx], df = df,cell_types = cell_types,seed = params$seed[i],add_count = params$add_count[i],phi_j = params$overdispersion[i],doublet_mode = T, my_beta=my_beta, X2=X2)
  res <- estimate_single_gene(myRCTD, params$gene[i],barcodes = data$barcodes,y = data$y,N = data$N,df = df,doublet_mode = T, my_beta=my_beta)
  metrics <- compute_metrics(res, data, est_col=ct_idx)
  cat(paste(c(i,round(metrics,4)), collapse=','), file = 'simulation_results_cere4.csv', sep='\n', append=T)
}
