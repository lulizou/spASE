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

PATH <- '../inst/extdata/slideseq/hippo3'
positions1 <- fread(file.path(PATH, 'Puck_221014_21_barcode_matching.txt.gz'), header=F, col.names = c('bead','idk1','x','y')) |>
  distinct(bead, .keep_all=T) |>
  mutate(bead = paste0(bead, '_1'))
counts1 <- fread(file.path(PATH, '2022-12-11_Puck_221014_21_combined_atropos_polyAGtrim_nm3.csv')) |>
  mutate(bead = paste0(substr(bead, 1, 14), '_1')) |>
  filter(bead %in% positions1$bead)
positions2 <- fread(file.path(PATH, 'Puck_221014_22_barcode_matching.txt.gz'), header=F, col.names = c('bead','idk1','x','y')) |>
  distinct(bead, .keep_all=T) |>
  mutate(bead = paste0(bead, '_2'))
counts2 <- fread(file.path(PATH, '2022-12-11_Puck_221014_22_combined_atropos_polyAGtrim_nm3.csv')) |>
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
# flip puck2 coordinates over x=2500
vline <- 2500
positions2 <- positions2 |> mutate(x = -(x-vline)+vline)
# shrink x coordinates
positions2 <- positions2 |> mutate(x = x/1.1)
# rotation
xcenter <- (max(positions2$x)+min(positions2$x))/2
ycenter <- (max(positions2$y)+min(positions2$y))/2
xshift <- 552
yshift <- 278
angle <- 47

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

saveRDS(myhippo, file = 'rctd_combined_hippo_3_nm3.rds')
