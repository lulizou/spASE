library(spacexr)
library(spASE)

res <- readRDS('results/rctd_mix_5_visium.rds')
maternal_counts_matrix <- res@originalSpatialRNA@maternalCounts
paternal_counts_matrix <- res@originalSpatialRNA@paternalCounts
myfit <- scase(
  maternal_counts_matrix,
  paternal_counts_matrix,
  min.cells = 2^7, cores=8, verbose=T
)
myfit <- scase(maternal_counts_matrix, paternal_counts_matrix, cores=1, verbose=T)
saveRDS(myfit, 'results/results_overall_bias_mix_5_visium.rds')
