library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(spASE)

PATH <- '../inst/extdata/slideseq/hippo1'
positions1 <- fread(file.path(PATH, 'Puck_200102_01_barcode_matching.txt'), col.names = c('idk1','idk2','bead','idk3','x','y')) |>
  distinct(bead, .keep_all=T) |>
  mutate(bead = paste0(bead, '_1'))
counts1 <- fread(file.path(PATH, '2020-08-20_Puck_200102_01_combined_aligned_counts_4-22-22.csv')) |>
  mutate(bead = paste0(substr(bead, 1, 14), '_1')) |>
  filter(bead %in% positions1$bead)
positions2 <- fread(file.path(PATH, 'Puck_200102_02_barcode_matching.txt'), col.names = c('idk1','idk2','bead','idk3','x','y')) |>
  distinct(bead, .keep_all=T) |>
  mutate(bead = paste0(bead, '_2'))
counts2 <- fread(file.path(PATH, '2020-08-20_Puck_200102_02_combined_aligned_counts_4-22-22.csv')) |>
  mutate(bead = paste0(substr(bead, 1, 14), '_2')) |>
  filter(bead %in% positions2$bead)

# MANUAL ALIGNING THE PUCKS
# flip puck2 coordinates over x=2500
angle <- 261
xshift <- 370
yshift <- -160
xcenter <- 3250
ycenter <- 3025

x <- positions2$x
y <- positions2$y
xnew <- (x-xcenter)*cos(angle*pi/180) - (y-ycenter)*sin(angle*pi/180)+xcenter + xshift
ynew <- (x-xcenter)*sin(angle*pi/180) + (y-ycenter)*cos(angle*pi/180)+ycenter + yshift

positions2$x <- xnew
positions2$y <- ynew

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


positions <- positions1 |> bind_rows(positions2) |>
  select(bead,x,y)

coords <- positions |> select(bead, x, y)

myfit <- spase(maternal_counts_matrix, paternal_counts_matrix, coords,  min.cells = 2^7,cores=8,verbose=T)
saveRDS(myfit, file = 'spasev1_results_global_combined_hippo_1.rds')

# add in cell type

rctd_results <- readRDS('rctd_combined_hippo_1.rds')
cell_types <- data.frame(bead = rownames(rctd_results@results$results_df), cell_type = rctd_results@results$results_df$first_type)
coords_cell_types <- coords |> left_join(cell_types, by = 'bead') |> filter(!is.na(cell_type))

myfit <- spase(maternal_counts_matrix, paternal_counts_matrix, coords_cell_types,cores=1,verbose=T)
saveRDS(myfit, file = 'spasev1_results_celltype_combined_hippo_1.rds')
