library(spacexr)
library(spASE)
#

hippo3 <- readRDS('rctd_combined_hippo_3_nm3.rds')
hippo3@config$max_cores <- 8
cell_types <- c('Astrocyte', 'CA1', 'CA3', 'Dentate', 'Interneuron', 'Oligodendrocyte')
hippo3_intercept <- run.CSIDE.intercept(hippo3, cell_types = cell_types, cell_type_threshold = 50, spase=F, logs=T)
hippo3_intercept <- run.CSIDE.intercept(hippo3_intercept, cell_types = cell_types, cell_type_threshold = 50, spase=T, logs=T)
save(hippo3, hippo3_intercept, file='cside_spase_combined_hippo_3_df_5.RData')
hippo3 <- run.CSIDE.nonparam(hippo3, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = F, logs=T,
                              test_genes_sig = T)
hippo3 <- run.CSIDE.nonparam(hippo3, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = T, logs=T,
                              test_genes_sig = T)
saveRDS(hippo3, file='cside_spase_combined_hippo_3_df_5.rds')
