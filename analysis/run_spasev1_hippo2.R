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
saveRDS(myfit, file = 'spasev1_results_global_hippo_2.rds')

# add in cell type

rctd_results <- readRDS('hippo_slideseq_mouse_2.rds')
res <- rctd_results@results$results_df
cell_types <- data.frame(bead = rownames(res), cell_type = res$first_type, spot_class = res$spot_class)
cell_types <- cell_types |> filter(spot_class != 'reject') |> select(bead, cell_type)
coords_cell_types <- coords |> left_join(cell_types, by = 'bead') |> filter(!is.na(cell_type))

myfit_celltype <- spase(maternal_counts1_matrix, paternal_counts1_matrix, coords_cell_types,min.pixels=2^7,cores=1,verbose=T)
saveRDS(myfit_celltype, file = 'spasev1_results_celltype_hippo_2.rds')
