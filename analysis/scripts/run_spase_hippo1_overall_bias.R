library(spacexr)
library(spASE)

res <- readRDS('results/rctd_hippo_1.rds')
maternal_counts_matrix <- res@originalSpatialRNA@maternalCounts
paternal_counts_matrix <- res@originalSpatialRNA@paternalCounts
myfit <- scase(
  maternal_counts_matrix,
  paternal_counts_matrix,
  min.cells = 2^7, cores=8, verbose=T
)
saveRDS(myfit, 'results_overall_bias_hippo_1.rds')
