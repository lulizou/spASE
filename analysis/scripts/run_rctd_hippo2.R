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

counts <- counts |>
  arrange(bead, gene) |>
  mutate(total = CAST+`129`+Unassigned) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
counts_matrix <- sparseMatrix(i=counts$gene, j=counts$bead, x=counts$total)
rownames(counts_matrix) <- levels(counts$gene)
colnames(counts_matrix) <- levels(counts$bead)


positions <- positions |>
  select(bead,x,y) |>
  column_to_rownames('bead')
puck <- SpatialRNA(positions, counts_matrix, maternalCounts = maternal_counts1_matrix, paternalCounts = paternal_counts1_matrix)
myhippo <- create.RCTD(puck, hippo_ref, max_cores=8)
myhippo <- run.RCTD(myhippo, doublet_mode = 'doublet')

saveRDS(myhippo, 'results/rctd_hippo_2.rds')

