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
saveRDS(myfit, 'results_overall_bias_hippo_2.rds')
