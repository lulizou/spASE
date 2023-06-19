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

maternal_counts1 <- counts |>
  arrange(bead, gene) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
maternal_counts1_matrix <- sparseMatrix(i=maternal_counts1$gene, j=maternal_counts1$bead, x=maternal_counts1$CAST)
rownames(maternal_counts1_matrix) <- levels(maternal_counts1$gene)
colnames(maternal_counts1_matrix) <- levels(maternal_counts1$bead)
paternal_counts1 <- counts |>
  arrange(bead, gene) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
paternal_counts1_matrix <- sparseMatrix(i=paternal_counts1$gene, j=paternal_counts1$bead, x=paternal_counts1$`129`)
rownames(paternal_counts1_matrix) <- levels(paternal_counts1$gene)
colnames(paternal_counts1_matrix) <- levels(paternal_counts1$bead)

positions <- positions |>
  select(barcode, starts_with('pxl')) |>
  dplyr::rename(x = pxl_row_in_fullres, y = pxl_col_in_fullres)


coords <- positions |> select(barcode, x, y)
myfit <- spase(maternal_counts1_matrix, paternal_counts1_matrix, coords,cores=8)
saveRDS(myfit, file = 'spase_results_visium_cere_3.rds')

# add in cell type

cere_vis <- readRDS('cerebellum_visium_mouse_3.rds')
results <- cere_vis@results
# normalize the cell type proportions to sum to 1.
norm_weights <- normalize_weights(results$weights)
cell_type_names <- cere_vis@cell_type_info$info[[2]] #list of cell type names
labels <- cell_type_names[max.col(norm_weights, 'first')]

cell_types <- data.frame(barcode = rownames(norm_weights), cell_type = labels)
coords_cell_types <- coords |> left_join(cell_types, by = 'barcode') |> filter(!is.na(cell_type))

myfit <- spase(maternal_counts1_matrix, paternal_counts1_matrix, coords_cell_types,cores=8,verbose=T)
saveRDS(myfit, file = 'spasev1_results_celltypecere_3.rds')
