library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)

PATH <- '../inst/extdata/slideseq/hippo3'
positions1 <- fread(file.path(PATH, 'Puck_221014_21_barcode_matching.txt.gz'), header=F, col.names = c('bead','idk1','x','y')) |>
  distinct(bead, .keep_all=T)
counts1 <- fread(file.path(PATH, '2022-12-11_Puck_221014_21_combined_atropos_polyAGtrim.csv')) |>
  mutate(bead = substr(bead, 1, 14)) |>
  filter(bead %in% positions1$bead)
positions2 <- fread(file.path(PATH, 'Puck_221014_22_barcode_matching.txt.gz'), header=F, col.names = c('bead','idk1','x','y')) |>
  distinct(bead, .keep_all=T)
counts2 <- fread(file.path(PATH, '2022-12-11_Puck_221014_22_combined_atropos_polyAGtrim.csv')) |>
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


myRCTD <- readRDS('hippo_slideseq_mouse_3_1.rds')
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
saveRDS(myRCTD, file.path('cside_spase_hippo_3_1.rds'))

myRCTD <- readRDS('hippo_slideseq_mouse_3_2.rds')
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
saveRDS(myRCTD, file.path('cside_spase_hippo_3_2.rds'))
