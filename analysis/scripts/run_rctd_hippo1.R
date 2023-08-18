library(spacexr)
library(spASE)
library(data.table)
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

counts <- counts1 |>
  bind_rows(counts2) |>
  arrange(bead, gene) |>
  mutate(total = CAST+`129`+Unassigned) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
counts_matrix <- sparseMatrix(i=counts$gene, j=counts$bead, x=counts$total)
rownames(counts_matrix) <- levels(counts$gene)
colnames(counts_matrix) <- levels(counts$bead)

# MANUAL ALIGNING THE PUCKS
# rotation
xcenter <- 3250
ycenter <- 3025
xshift <- 370
yshift <- -160
angle <- 261*pi/180
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
  select(bead,x,y) |>
  column_to_rownames('bead')
counts <- rbind(counts1, counts2)

# First run RCTD

puck <- SpatialRNA(positions, counts_matrix, maternalCounts = maternal_counts_matrix, paternalCounts = paternal_counts_matrix)
myhippo <- create.RCTD(puck, hippo_ref, max_cores=8)
myhippo <- run.RCTD(myhippo, doublet_mode = 'doublet')

saveRDS(myhippo, file = 'results/rctd_hippo_1.rds')
