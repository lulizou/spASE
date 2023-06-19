library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)

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


myRCTD <- readRDS('hippo_slideseq_mouse_1_1.rds')
myRCTD@config$max_cores <- 8
myRCTD@spatialRNA@maternalCounts <- maternal_counts1_matrix
myRCTD@spatialRNA@paternalCounts <- paternal_counts1_matrix
myRCTD@originalSpatialRNA@maternalCounts <- maternal_counts1_matrix
myRCTD@originalSpatialRNA@paternalCounts <- paternal_counts1_matrix
myRCTD@spase_results <- list()
cell_types <- c('Astrocyte', 'CA1', 'CA3', 'Dentate', 'Interneuron', 'Oligodendrocyte')
myRCTD <- run.CSIDE.nonparam(myRCTD, df = 5, cell_types = cell_types,
                             cell_type_threshold = 50, spase = F)
myRCTD <- run.CSIDE.nonparam(myRCTD, df = 5, cell_types = cell_types,
                             cell_type_threshold = 50, spase = T)
saveRDS(myRCTD, file.path('cside_spase_hippo_1_1.rds'))

myRCTD <- readRDS('hippo_slideseq_mouse_1_2.rds')
myRCTD@config$max_cores <- 8
myRCTD@spatialRNA@maternalCounts <- maternal_counts2_matrix
myRCTD@spatialRNA@paternalCounts <- paternal_counts2_matrix
myRCTD@originalSpatialRNA@maternalCounts <- maternal_counts2_matrix
myRCTD@originalSpatialRNA@paternalCounts <- paternal_counts2_matrix
myRCTD@spase_results <- list()
cell_types <- c('Astrocyte', 'CA1', 'CA3', 'Dentate', 'Interneuron', 'Oligodendrocyte')
myRCTD <- run.CSIDE.nonparam(myRCTD, df = 5, cell_types = cell_types,
                             cell_type_threshold = 50, spase = F)
myRCTD <- run.CSIDE.nonparam(myRCTD, df = 5, cell_types = cell_types,
                             cell_type_threshold = 50, spase = T)
saveRDS(myRCTD, file.path('cside_spase_hippo_1_2.rds'))
