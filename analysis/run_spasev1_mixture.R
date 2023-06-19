library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(spASE)


PATH <- '../inst/extdata/visium/mixture'
counts <- fread(file.path(PATH, 'mixture.csv'))
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
saveRDS(myfit, file = 'spase_results_visium_mixture_5.rds')
