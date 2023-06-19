library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(spASE)

PATH <- '../inst/extdata/slideseq/cere4'

positions1 <- fread(file.path(PATH, 'Puck_221014_23_barcode_matching.txt.gz'), header=F, col.names = c('bead','idk1','x','y')) |>
  distinct(bead, .keep_all=T) |>
  mutate(bead = paste0(bead, '_1'))
counts1 <- fread(file.path(PATH, '2022-12-11_Puck_221014_23_combined_atropos_polyAGtrim_nm3.csv')) |>
  mutate(bead = paste0(substr(bead, 1, 14), '_1')) |>
  filter(bead %in% positions1$bead)
positions2 <- fread(file.path(PATH, 'Puck_221014_24_barcode_matching.txt.gz'), header=F, col.names = c('bead','idk1','x','y')) |>
  distinct(bead, .keep_all=T) |>
  mutate(bead = paste0(bead, '_2'))
counts2 <- fread(file.path(PATH, '2022-12-11_Puck_221014_24_combined_atropos_polyAGtrim_nm3.csv')) |>
  mutate(bead = paste0(substr(bead, 1, 14), '_2')) |>
  filter(bead %in% positions2$bead)

# MANUAL ALIGNING THE PUCKS
# rotation
xcenter <- (max(positions2$x)+min(positions2$x))/2
ycenter <- (max(positions2$y)+min(positions2$y))/2
xshift <- 748
yshift <- -241
angle <- 86

x <- positions2$x
y <- positions2$y
xnew <- (x-xcenter)*cos(angle*pi/180) - (y-ycenter)*sin(angle*pi/180)+xcenter + xshift
ynew <- (x-xcenter)*sin(angle*pi/180) + (y-ycenter)*cos(angle*pi/180)+ycenter + yshift

positions2$x <- xnew
positions2$y <- ynew

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

myfit <- spase(maternal_counts_matrix, paternal_counts_matrix, coords,cores=8,verbose=T)
saveRDS(myfit, file = 'spasev1_results_global_combined_cere_4_nm3.rds')

# add in cell type

# rctd_results <- readRDS('rctd_combined_cere_4.rds')
# res <- rctd_results@results$results_df
# cell_types <- data.frame(bead = rownames(res), cell_type = res$first_type, spot_class = res$spot_class)
# cell_types <- cell_types |> filter(spot_class != 'reject') |> select(bead, cell_type)
# coords_cell_types <- coords |> left_join(cell_types, by = 'bead') |> filter(!is.na(cell_type))
#
# myfit <- spase(maternal_counts_matrix, paternal_counts_matrix, coords_cell_types,min.pixels=2^7,cores=8,verbose=T)
# saveRDS(myfit, file = 'spasev1_results_celltype_combined_cere_4.rds')
