library(spacexr)
library(spASE)


cere3 <- readRDS('results/rctd_cere_3.rds')
cere3@config$max_cores <- 8
cell_types <- c('Fibroblast', 'Granule', 'MLI2', 'Purkinje', 'Bergmann', 'Oligodendrocytes')
cere3_intercept <- run.CSIDE.intercept(cere3, cell_types = cell_types, cell_type_threshold = 50, spase=F, logs=T)
cere3_intercept <- run.CSIDE.intercept(cere3_intercept, cell_types = cell_types, cell_type_threshold = 50, spase=T, logs=T)
saveRDS(cere3_intercept, file='results/results_celltype_intercept_cere_3_df_5.rds')
cere3 <- run.CSIDE.nonparam(cere3, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = F, logs=T,
                             test_genes_sig = T)
cere3 <- run.CSIDE.nonparam(cere3, df = 5, cell_types = cell_types,
                              cell_type_threshold = 2^7, spase = T, logs=T,
                             test_genes_sig = T)
saveRDS(cere3, file='results/results_celltype_cere_3_df_5.rds')
