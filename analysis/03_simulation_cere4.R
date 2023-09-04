# first argument is the chunk of params to run (1-18)
library(spacexr)
library(spASE)
devtools::load_all()

args <- commandArgs(trailingOnly=T)

source('../R/spase_sim_utils.R')


myRCTD <- readRDS('results/results_celltype_cere_4_visium_df_5.rds')
params <- read.delim('simulations_cere4_params.tsv')
df <- 5
cell_types <- colnames(myRCTD@de_results$gene_fits$mean_val)

chunk <- as.integer(args[1])
idx <- (1:1000)+(chunk-1)*1000
print(args)
print(idx)
for (i in idx) {
  message((paste(Sys.time(), i)))
  ct_idx <- which(cell_types == params$cell_type[i])
  data <- simulate_data(myRCTD, params$gene[i], cell_types[ct_idx], df = df,cell_types = cell_types,seed = params$seed[i],add_count = params$add_count[i],phi_j = params$overdispersion[i],doublet_mode = F)
  res <- estimate_single_gene(myRCTD, params$gene[i],barcodes = data$barcodes,y = data$y,N = data$N,df = df,doublet_mode = F)
  metrics <- compute_metrics(res, data, est_col=ct_idx)
  cat(paste(c(i,round(metrics,4)), collapse=','), file = 'simulation_results_cere4.csv', sep='\n', append=T)
}
