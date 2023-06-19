library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(spASE)

PATH <- '../inst/extdata/slideseq/cere1'

positions1 <- fread(file.path(PATH, 'Puck_221014_23_barcode_matching.txt.gz'), header=F, col.names = c('bead','idk1','x','y')) |>
  distinct(bead, .keep_all=T)
counts1 <- fread(file.path(PATH, '2022-12-11_Puck_221014_23_combined_atropos_polyAGtrim.csv')) |>
  mutate(bead = substr(bead, 1, 14)) |>
  filter(bead %in% positions1$bead)
positions2 <- fread(file.path(PATH, 'Puck_221014_24_barcode_matching.txt.gz'), header=F, col.names = c('bead','idk1','x','y')) |>
  distinct(bead, .keep_all=T)
counts2 <- fread(file.path(PATH, '2022-12-11_Puck_221014_24_combined_atropos_polyAGtrim.csv')) |>
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
myfit <- spase(maternal_counts1_matrix, paternal_counts1_matrix, coords,cores=8,verbose=T)
saveRDS(myfit, file = 'spase_results_cere_4_1.rds')

coords <- positions2 |> select(bead, x, y)
myfit <- spase(maternal_counts2_matrix, paternal_counts2_matrix, coords,cores=8,verbose=T)
saveRDS(myfit, file = 'spase_results_cere_4_2.rds')


