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


myfit <- scase(maternal_counts_matrix, paternal_counts_matrix, cores=1, verbose=T)
saveRDS(myfit, 'betabin_results_overall_bias_mixture_visium.rds')
