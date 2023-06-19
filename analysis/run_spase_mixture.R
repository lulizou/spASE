library(data.table)
library(spacexr)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(spASE)


PATH <- '../inst/extdata/visium/mixture'
counts <- fread(file.path(PATH, 'mixture.csv'))
positions <- fread(file.path(PATH, 'tissue_positions.csv'))

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

myRCTD <- readRDS('mixture_visium_cerebellumRef_mouse_5.rds')
myRCTD@config$max_cores <- 8
myRCTD@spatialRNA@maternalCounts <- maternal_counts1_matrix
myRCTD@spatialRNA@paternalCounts <- paternal_counts1_matrix
myRCTD@originalSpatialRNA@maternalCounts <- maternal_counts1_matrix
myRCTD@originalSpatialRNA@paternalCounts <- paternal_counts1_matrix
myRCTD@spase_results <- list()
myRCTD@config$doublet_mode <- 'full'
cell_types <- c('Fibroblast', 'Granule', 'MLI2', 'Purkinje', 'Bergmann', 'Oligodendrocytes')
myRCTD <- run.CSIDE.nonparam(myRCTD, df = 5, cell_types = cell_types,
                             cell_type_threshold = 50, spase = F, doublet_mode = F)
myRCTD <- run.CSIDE.nonparam(myRCTD, df = 5, cell_types = cell_types,
                             cell_type_threshold = 50, spase = T, doublet_mode = F)
saveRDS(myRCTD, file.path('cside_spase_mixture_5.rds'))
