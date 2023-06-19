
library(spacexr)
library(spASE)

hippo2 <- readRDS('hippo_slideseq_mouse_2.rds')
hippo2@config$max_cores <- 8
cell_types <- c('Astrocyte', 'CA1', 'CA3', 'Dentate', 'Interneuron', 'Oligodendrocyte')
hippo2_intercept <- run.CSIDE.intercept(hippo2, cell_types = cell_types, cell_type_threshold = 50, spase=F, logs=T)
hippo2_intercept <- run.CSIDE.intercept(hippo2_intercept, cell_types = cell_types, cell_type_threshold = 50, spase=T, logs=T)
save(hippo2, hippo2_intercept, file='cside_spase_combined_hippo_2_df_5.RData')
hippo2 <- run.CSIDE.nonparam(hippo2, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = F, logs=T,
                             test_genes_sig = T)
hippo2 <- run.CSIDE.nonparam(hippo2, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = T, logs=T,
                             test_genes_sig = T)
saveRDS(hippo2, file='cside_spase_combined_hippo_2_df_5.rds')
