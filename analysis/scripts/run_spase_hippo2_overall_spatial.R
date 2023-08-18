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

coords <- positions |> select(bead, x, y)
myfit <- spase(maternal_counts1_matrix, paternal_counts1_matrix, coords, min.pixels=2^7, cores=8)
saveRDS(myfit, file = 'results_overall_spatial_hippo_2.rds')
