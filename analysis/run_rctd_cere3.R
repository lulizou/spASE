library(spacexr)
library(spASE)
library(tibble)


cerebellum_ref <- readRDS('../inst/extdata/reference_scrna/cerebellum_sn_ref_spacexr.rds')
PATH <- '../inst/extdata/visium/cerebellum_spatial'
counts <- fread(file.path(PATH, 'cerebellum.csv'))
positions <- fread(file.path(PATH, 'tissue_positions.csv'))


counts <- counts |>
  arrange(bead, gene) |>
  mutate(total = CAST+`129`+Unassigned) |>
  mutate(bead = as.factor(bead), gene = as.factor(gene))
counts_matrix <- sparseMatrix(i=counts$gene, j=counts$bead, x=counts$total)
rownames(counts_matrix) <- levels(counts$gene)
colnames(counts_matrix) <- levels(counts$bead)


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

positions <- positions |>
  select(barcode, starts_with('pxl')) |>
  column_to_rownames('barcode') |>
  rename(x = pxl_row_in_fullres, y = pxl_col_in_fullres)
puck <- SpatialRNA(positions, counts_matrix, maternalCounts = maternal_counts1_matrix,
                   paternalCounts = paternal_counts1_matrix)
myCerebellum <- create.RCTD(puck, cerebellum_ref, max_cores=8)
myCerebellum <- run.RCTD(myCerebellum, doublet_mode = 'full')

saveRDS(myCerebellum, 'rctd_cere_3.rds')
