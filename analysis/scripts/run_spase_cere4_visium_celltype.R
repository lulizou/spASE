library(spacexr)
library(spASE)


cere4 <- readRDS('rctd_cere_4_visium.rds')
cere4@config$max_cores <- 8
cell_types <- c('Fibroblast', 'Granule', 'MLI2','Oligodendrocytes')
cere4_intercept <- run.CSIDE.intercept(cere4, cell_types = cell_types,
                                        cell_type_threshold = 2^7, spase=F,
                                        logs=T, doublet_mode = F)
 cere4_intercept <- run.CSIDE.intercept(cere4_intercept, cell_types = cell_types,
                                        cell_type_threshold = 2^7, spase=T, logs=T, doublet_mode = F)
 saveRDS(cere4_intercept, file='results_celltype_intercept_cere_4_visium_df_5.rds')
cere4 <- run.CSIDE.nonparam(cere4, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = F, logs=T,
                             test_genes_sig = T, doublet_mode = F)
cere4 <- run.CSIDE.nonparam(cere4, df = 5, cell_types = cell_types,
                             cell_type_threshold = 2^7, spase = T, logs=T,
                             test_genes_sig = T, doublet_mode = F)
saveRDS(cere4, file='results_celltype_cere_4_visium_df_5.rds')
