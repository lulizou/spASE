

library(spacexr)
library(spASE)

hippo1 <- readRDS('rctd_combined_hippo_1.rds')
hippo1@config$max_cores <- 8
cell_types <- c('Astrocyte', 'CA1', 'CA3', 'Dentate', 'Interneuron', 'Oligodendrocyte')
hippo1_intercept <- run.CSIDE.intercept(hippo1, cell_types = cell_types, cell_type_threshold = 50, spase=F, logs=T)
hippo1_intercept <- run.CSIDE.intercept(hippo1_intercept, cell_types = cell_types, cell_type_threshold = 50, spase=T, logs=T)
save(hippo1, hippo1_intercept, file='cside_spase_combined_hippo_1_df_5.RData')
hippo1 <- run.CSIDE.nonparam(hippo1, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = F, logs=T,
                              test_genes_sig = T)
hippo1 <- run.CSIDE.nonparam(hippo1, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = T, logs=T,
                              test_genes_sig = T)
saveRDS(hippo1, file='cside_spase_combined_hippo_1_df_5.rds')
