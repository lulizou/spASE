library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(spASE)

PATH <- '../inst/extdata/slideseq/hippo2'

positions <- fread(file.path(PATH, 'Puck_201006_22_barcode_matching.txt'), col.names = c('idk1','idk2','bead','idk3','x','y'))  |>
  distinct(bead, .keep_all=T)
counts <- fread(file.path(PATH, '2020-11-25_Puck_201006_22_combined_aligned_k100_local.csv')) |>
  mutate(bead = substr(bead, 1, 14)) |>
  filter(bead %in% positions$bead)

maternal_counts <- counts |>
  arrange(bead, gene) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
maternal_counts_matrix <- sparseMatrix(i=maternal_counts$gene, j=maternal_counts$bead, x=maternal_counts$CAST)
rownames(maternal_counts_matrix) <- levels(maternal_counts$gene)
colnames(maternal_counts_matrix) <- levels(maternal_counts$bead)
paternal_counts <- counts |>
  arrange(bead, gene) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
paternal_counts_matrix <- sparseMatrix(i=paternal_counts$gene, j=paternal_counts$bead, x=paternal_counts$`129`)
rownames(paternal_counts_matrix) <- levels(paternal_counts$gene)
colnames(paternal_counts_matrix) <- levels(paternal_counts$bead)

myfit <- scase(maternal_counts_matrix, paternal_counts_matrix, min.cells = 2^7, cores=8, verbose=T)
saveRDS(myfit, 'betabin_results_overall_bias_hippo_2.rds')

rctd_results <- readRDS('hippo_slideseq_mouse_2.rds')
res <- rctd_results@results$results_df
cell_types <- data.frame(bead = rownames(res), cell_type = res$first_type, spot_class = res$spot_class)
cell_types <- cell_types |> filter(spot_class != 'reject') |> select(bead, cell_type)

maternal_counts_bead <- maternal_counts_matrix[,cell_types$bead]
paternal_counts_bead <- paternal_counts_matrix[,cell_types$bead]

myfit_celltype <- scase(maternal_counts_bead, paternal_counts_bead,
                        covariates = cell_types, min.cells = 2^7, cores=8, verbose=T)
saveRDS(myfit_celltype, 'betabin_results_celltype_bias_hippo_2.rds')
