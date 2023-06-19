library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(spASE)

PATH <- '../inst/extdata/slideseq/cere4'

positions1 <- fread(file.path(PATH, 'Puck_221014_23_barcode_matching.txt.gz'), header=F, col.names = c('bead','idk1','x','y')) |>
  distinct(bead, .keep_all=T) |>
  mutate(bead = paste0(bead, '_1'))
counts1 <- fread(file.path(PATH, '2022-12-11_Puck_221014_23_combined_atropos_polyAGtrim_nm3.csv')) |>
  mutate(bead = paste0(substr(bead, 1, 14), '_1')) |>
  filter(bead %in% positions1$bead)
positions2 <- fread(file.path(PATH, 'Puck_221014_24_barcode_matching.txt.gz'), header=F, col.names = c('bead','idk1','x','y')) |>
  distinct(bead, .keep_all=T) |>
  mutate(bead = paste0(bead, '_2'))
counts2 <- fread(file.path(PATH, '2022-12-11_Puck_221014_24_combined_atropos_polyAGtrim_nm3.csv')) |>
  mutate(bead = paste0(substr(bead, 1, 14), '_2')) |>
  filter(bead %in% positions2$bead)

maternal_counts<- counts1 |>
  bind_rows(counts2) |>
  arrange(bead, gene) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
maternal_counts_matrix <- sparseMatrix(i=maternal_counts$gene, j=maternal_counts$bead, x=maternal_counts$CAST)
rownames(maternal_counts_matrix) <- levels(maternal_counts$gene)
colnames(maternal_counts_matrix) <- levels(maternal_counts$bead)

paternal_counts <- counts1 |>
  bind_rows(counts2) |>
  arrange(bead, gene) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
paternal_counts_matrix <- sparseMatrix(i=paternal_counts$gene, j=paternal_counts$bead, x=paternal_counts$`129`)
rownames(paternal_counts_matrix) <- levels(paternal_counts$gene)
colnames(paternal_counts_matrix) <- levels(paternal_counts$bead)

myfit <- scase(maternal_counts_matrix, paternal_counts_matrix, min.cells = 2^7, cores=8, verbose=T)
saveRDS(myfit, 'betabin_results_overall_bias_cere_4_nm3.rds')

rctd_results <- readRDS('rctd_combined_cere_4.rds')
res <- rctd_results@results$results_df
cell_types <- data.frame(bead = rownames(res), cell_type = res$first_type, spot_class = res$spot_class)
cell_types <- cell_types |> filter(spot_class != 'reject') |> select(bead, cell_type)

maternal_counts_bead <- maternal_counts_matrix[,cell_types$bead]
paternal_counts_bead <- paternal_counts_matrix[,cell_types$bead]

myfit_celltype <- scase(maternal_counts_bead, paternal_counts_bead,
                        covariates = cell_types, min.cells = 2^7, cores=8, verbose=T)
saveRDS(myfit_celltype, 'betabin_results_celltype_bias_cere_4_nm3.rds')
