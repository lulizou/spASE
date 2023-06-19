library(spacexr)
library(spASE)


cere3 <- readRDS('rctd_cere_3.rds')
cere3@config$max_cores <- 8
cell_types <- c('Fibroblast', 'Granule', 'MLI2','Oligodendrocytes')
cere3_intercept <- run.CSIDE.intercept(cere3, cell_types = cell_types,
                                        cell_type_threshold = 2^7, spase=F,
                                        logs=T, doublet_mode = F)
 cere3_intercept <- run.CSIDE.intercept(cere3_intercept, cell_types = cell_types,
                                        cell_type_threshold = 2^7, spase=T, logs=T, doublet_mode = F)
 save(cere3, cere3_intercept, file='cside_spase_combined_cere_3_df_5.RData')
 cere3 <- run.CSIDE.nonparam(cere3, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = F, logs=T,
                             test_genes_sig = T, doublet_mode = F)
 saveRDS(cere3, file='cside_spase_combined_cere_3_df_5.rds')
cell_types <- c('Fibroblast', 'Granule', 'MLI2','Oligodendrocytes')
cere3 <- run.CSIDE.nonparam(cere3, df = 5, cell_types = cell_types,
                             cell_type_threshold = 2^7, spase = T, logs=T,
                             test_genes_sig = T, doublet_mode = F)
saveRDS(cere3, file='cside_spase_combined_cere_3_df_5.rds')
