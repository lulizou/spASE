library(spacexr)
library(spASE)


cere4 <- readRDS('rctd_combined_cere_4_nm3.rds')
cere4@config$max_cores <- 8
cell_types <- c('Fibroblast', 'Granule', 'MLI2', 'Purkinje', 'Bergmann', 'Oligodendrocytes')
cere4_intercept <- run.CSIDE.intercept(cere4, cell_types = cell_types, cell_type_threshold = 50, spase=F, logs=T)
cere4_intercept <- run.CSIDE.intercept(cere4_intercept, cell_types = cell_types, cell_type_threshold = 50, spase=T, logs=T)
save(cere4, cere4_intercept, file='cside_spase_combined_cere_4_df_5.RData')
cere4 <- run.CSIDE.nonparam(cere4, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = F, logs=T,
                             test_genes_sig = T)
saveRDS(cere4, file='cside_spase_combined_cere_4_df_5.rds')
cere4 <- run.CSIDE.nonparam(cere4, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = T, logs=T,
                             test_genes_sig = T)
saveRDS(cere4, file='cside_spase_combined_cere_4_df_5.rds')
