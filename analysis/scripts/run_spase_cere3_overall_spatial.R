library(spacexr)
library(spASE)

res <- readRDS('results/rctd_cere_3.rds')
maternal_counts_matrix <- res@originalSpatialRNA@maternalCounts
paternal_counts_matrix <- res@originalSpatialRNA@paternalCounts
coords <- res@originalSpatialRNA@coords
coords$bead <- rownames(coords)
coords <- coords[c('bead','x','y')]

myfit <- spase(maternal_counts_matrix, paternal_counts_matrix, coords,cores=8,verbose=T)
saveRDS(myfit, file = 'results/results_overall_spatial_cere_3.rds')


