library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(spASE)

PATH <- '../inst/extdata/visium/cerebellum_spatial'
counts <- fread(file.path(PATH, 'cerebellum.csv'))
positions <- fread(file.path(PATH, 'tissue_positions.csv'))

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

positions <- positions |>
  select(barcode, starts_with('pxl')) |>
  dplyr::rename(x = pxl_row_in_fullres, y = pxl_col_in_fullres)


coords <- positions |> select(barcode, x, y)


myfit <- scase(maternal_counts_matrix, paternal_counts_matrix, cores=8)
saveRDS(myfit, 'betabin_results_overall_bias_cere_visium.rds')

rctd_results <- readRDS('cerebellum_visium_mouse_3.rds')
results <- rctd_results@results
# normalize the cell type proportions to sum to 1.
norm_weights <- normalize_weights(results$weights)
cell_type_names <- rctd_results@cell_type_info$info[[2]] #list of cell type names
labels <- cell_type_names[max.col(norm_weights, 'first')]

cell_types <- data.frame(bead = rownames(norm_weights), cell_type = labels)
cell_types$cell_type <- factor(cell_types$cell_type)

maternal_counts_bead <- maternal_counts_matrix[,cell_types$bead]
paternal_counts_bead <- paternal_counts_matrix[,cell_types$bead]

myfit_celltype <- scase(maternal_counts_bead, paternal_counts_bead,
                        covariates = cell_types, cores=8, verbose=T)
saveRDS(myfit_celltype, 'betabin_results_celltype_bias_cere_visium.rds')
