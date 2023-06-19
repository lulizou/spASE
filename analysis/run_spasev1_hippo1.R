library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(spASE)

PATH <- '../inst/extdata/slideseq/hippo1'
positions1 <- fread(file.path(PATH, 'Puck_200102_01_barcode_matching.txt'), col.names = c('idk1','idk2','bead','idk3','x','y')) |>
  distinct(bead, .keep_all=T)
counts1 <- fread(file.path(PATH, '2020-08-20_Puck_200102_01_combined_aligned_counts_4-22-22.csv')) |>
  mutate(bead = substr(bead, 1, 14)) |>
  filter(bead %in% positions1$bead)
positions2 <- fread(file.path(PATH, 'Puck_200102_02_barcode_matching.txt'), col.names = c('idk1','idk2','bead','idk3','x','y')) |>
  mutate(bead = substr(bead, 1, 14)) |>
  distinct(bead, .keep_all=T)
counts2 <- fread(file.path(PATH, '2020-08-20_Puck_200102_02_combined_aligned_counts_4-22-22.csv')) |>
  mutate(bead = substr(bead, 1, 14)) |>
  filter(bead %in% positions2$bead)

maternal_counts1 <- counts1 |>
  arrange(bead, gene) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
maternal_counts1_matrix <- sparseMatrix(i=maternal_counts1$gene, j=maternal_counts1$bead, x=maternal_counts1$CAST)
rownames(maternal_counts1_matrix) <- levels(maternal_counts1$gene)
colnames(maternal_counts1_matrix) <- levels(maternal_counts1$bead)
paternal_counts1 <- counts1 |>
  arrange(bead, gene) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
paternal_counts1_matrix <- sparseMatrix(i=paternal_counts1$gene, j=paternal_counts1$bead, x=paternal_counts1$`129`)
rownames(paternal_counts1_matrix) <- levels(paternal_counts1$gene)
colnames(paternal_counts1_matrix) <- levels(paternal_counts1$bead)

maternal_counts2 <- counts2 |>
  arrange(bead, gene) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
maternal_counts2_matrix <- sparseMatrix(i=maternal_counts2$gene, j=maternal_counts2$bead, x=maternal_counts2$CAST)
rownames(maternal_counts2_matrix) <- levels(maternal_counts2$gene)
colnames(maternal_counts2_matrix) <- levels(maternal_counts2$bead)
paternal_counts2 <- counts2 |>
  arrange(bead, gene) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
paternal_counts2_matrix <- sparseMatrix(i=paternal_counts2$gene, j=paternal_counts2$bead, x=paternal_counts2$`129`)
rownames(paternal_counts2_matrix) <- levels(paternal_counts2$gene)
colnames(paternal_counts2_matrix) <- levels(paternal_counts2$bead)

coords <- positions1 |> select(bead, x, y)
myfit <- spase(maternal_counts1_matrix, paternal_counts1_matrix, coords,cores=8)
saveRDS(myfit, file = 'spase_results_hippo_1_1.rds')

coords <- positions2 |> select(bead, x, y)
myfit <- spase(maternal_counts2_matrix, paternal_counts2_matrix, coords,cores=8)
saveRDS(myfit, file = 'spase_results_hippo_1_2.rds')


