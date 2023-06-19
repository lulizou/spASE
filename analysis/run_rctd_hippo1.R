library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)

hippo_ref <- readRDS('../inst/extdata/reference_scrna/scRefSubsampled1000.rds')

celltypes_hippo <- hippo_ref@meta.data$liger_ident_coarse
levels(celltypes_hippo)[levels(celltypes_hippo)=='Denate'] <- 'Dentate'
hippo_counts <- hippo_ref@assays$RNA@counts
names(celltypes_hippo) <- colnames(hippo_counts)

hippo_ref <- Reference(counts = hippo_counts, cell_types = celltypes_hippo)

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


counts1 <- counts1 |>
  arrange(bead, gene) |>
  mutate(total = CAST+`129`+Unassigned) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
counts1_matrix <- sparseMatrix(i=counts1$gene, j=counts1$bead, x=counts1$total)
rownames(counts1_matrix) <- levels(counts1$gene)
colnames(counts1_matrix) <- levels(counts1$bead)

counts2 <- counts2 |>
  arrange(bead, gene) |>
  mutate(total = CAST+`129`+Unassigned) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
counts2_matrix <- sparseMatrix(i=counts2$gene, j=counts2$bead, x=counts2$total)
rownames(counts2_matrix) <- levels(counts2$gene)
colnames(counts2_matrix) <- levels(counts2$bead)

positions1 <- positions1 |>
  select(bead,x,y) |>
  column_to_rownames('bead')
positions2 <- positions2 |>
  select(bead,x,y) |>
  column_to_rownames('bead')
puck1 <- SpatialRNA(positions1, counts1_matrix)
myhippo1 <- create.RCTD(puck1, hippo_ref, max_cores=8)
myhippo1 <- run.RCTD(myhippo1, doublet_mode = 'doublet')

saveRDS(myhippo1, 'hippo_slideseq_mouse_1_1.rds')

puck2 <- SpatialRNA(positions2, counts2_matrix)
myhippo2 <- create.RCTD(puck2, hippo_ref, max_cores=8)
myhippo2 <- run.RCTD(myhippo2, doublet_mode = 'doublet')

saveRDS(myhippo2, 'hippo_slideseq_mouse_1_2.rds')
