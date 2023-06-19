library(dplyr)
library(tidyr)

gencode <- import('/rafalab/lzou/resources/gencode.vM10.annotation.gff3.gz')
xchr_genes <- unique(gencode$gene_name[which(seqnames(gencode)=='chrX')])
xchr_genes <- c(xchr_genes, 'Bex3')

# Visium

overall_bias_df <- readRDS('betabin_results_overall_bias_cere_visium.rds') |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))


# Total number of genes

overall_bias_df |>
  filter(totalUMI > 2^7) |>
  group_by(xchr) |>
  summarise(n())

visium_all_genes <- overall_bias_df |>
  filter(totalUMI > 2^7) |>
  select(gene,xchr)

# Overall p \neq 0.5

overall_bias_df |>
  filter(totalUMI > 2^7, qval < 0.01, p > 0.6 | p < 0.4) |>
  group_by(xchr) |>
  summarise(n())

visium_bias_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, qval < 0.01, p > 0.6 | p < 0.4) |>
  select(gene,xchr)

# Within cell type p \neq 0.5

celltype_bias_df <- readRDS('betabin_results_celltype_bias_cere_visium.rds') |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

celltype_bias_df |>
  filter(p > 0.6 | p < 0.4) |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

visium_ct_bias_genes <- celltype_bias_df |>
  filter(p > 0.6 | p < 0.4) |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)


# Overall spatial pattern

overall_spatial_df <- readRDS('spasev1_results_visium_cere_3.rds')$result |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

overall_spatial_df |>
  filter(qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

visium_overall_spatial_genes <- overall_spatial_df |>
  filter(qval < 0.01) |>
  select(gene,xchr)

# Within cell type spatial pattern

celltype_spatial_df <- readRDS('spasev1_results_celltypecere_3.rds')$result |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

celltype_spatial_df |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

visium_ct_spatial_genes <- celltype_spatial_df |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Slide-seq

overall_bias_df <- readRDS('betabin_results_overall_bias_cere_4.rds') |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))


# Total number of genes

overall_bias_df |>
  filter(totalUMI > 2^7) |>
  group_by(xchr) |>
  summarise(n())

slideseq_all_genes <- overall_bias_df |>
  filter(totalUMI > 2^7) |>
  select(gene,xchr)

# Overall p \neq 0.5

overall_bias_df |>
  filter(totalUMI > 2^7, qval < 0.01, p > 0.6 | p < 0.4) |>
  group_by(xchr) |>
  summarise(n())

slideseq_bias_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, qval < 0.01, p > 0.6 | p < 0.4) |>
  select(gene,xchr)

# Within cell type p \neq 0.5

celltype_bias_df <- readRDS('betabin_results_celltype_bias_cere_4.rds') |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

celltype_bias_df |>
  filter(p > 0.6 | p < 0.4) |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

slideseq_ct_bias_genes <- celltype_bias_df |>
  filter(p > 0.6 | p < 0.4) |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Overall spatial pattern

overall_spatial_df <- readRDS('spasev1_results_global_combined_cere_4.rds')$result |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

overall_spatial_df |>
  filter(qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

slideseq_overall_spatial_genes <- overall_spatial_df |>
  filter(qval < 0.01) |>
  select(gene,xchr)

# Within cell type spatial pattern

celltype_spatial_df <- readRDS('spasev1_results_celltype_combined_cere_4.rds')$result |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

celltype_spatial_df |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

slideseq_ct_spatial_genes <- celltype_spatial_df |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)


# Overlap

# Total genes

visium_all_genes |>
  mutate(in_slideseq = ifelse(gene %in% slideseq_all_genes$gene, 1, 0)) |>
  group_by(xchr) |>
  summarise(sum(in_slideseq))

# No significant ASE

visium_all_genes |>
  filter(!gene %in% c(visium_bias_genes$gene, visium_ct_bias_genes$gene,
                      visium_overall_spatial_genes$gene, visium_ct_spatial_genes$gene)) |>
  group_by(xchr) |>
  summarise(n())

slideseq_all_genes |>
  filter(!gene %in% c(slideseq_bias_genes$gene, slideseq_ct_bias_genes$gene,
                      slideseq_overall_spatial_genes$gene, slideseq_ct_spatial_genes$gene)) |>
  group_by(xchr) |>
  summarise(n())

visium_no_ase_genes <- visium_all_genes |>
  filter(!gene %in% c(visium_bias_genes$gene, visium_ct_bias_genes$gene,
                      visium_overall_spatial_genes$gene, visium_ct_spatial_genes$gene)) |>
  select(gene, xchr)

slideseq_no_ase_genes <- slideseq_all_genes |>
  filter(!gene %in% c(slideseq_bias_genes$gene, slideseq_ct_bias_genes$gene,
                      slideseq_overall_spatial_genes$gene, slideseq_ct_spatial_genes$gene)) |>
  select(gene, xchr)

visium_no_ase_genes |>
  mutate(in_slideseq = ifelse(gene %in% slideseq_no_ase_genes$gene, 1, 0)) |>
  group_by(xchr) |>
  summarise(sum(in_slideseq))

# Overall bias

visium_bias_genes |>
  mutate(in_slideseq = ifelse(gene %in% slideseq_bias_genes$gene, 1, 0)) |>
  group_by(xchr) |>
  summarise(sum(in_slideseq))

# Cell type bias

visium_ct_bias_genes |>
  mutate(in_slideseq = ifelse(gene %in% slideseq_ct_bias_genes$gene, 1, 0)) |>
  group_by(xchr) |>
  summarise(sum(in_slideseq))

# Overall spatial

visium_overall_spatial_genes |>
  mutate(in_slideseq = ifelse(gene %in% slideseq_overall_spatial_genes$gene, 1, 0)) |>
  group_by(xchr) |>
  summarise(sum(in_slideseq))

# Cell type spatial

visium_ct_spatial_genes |>
  mutate(in_slideseq = ifelse(gene %in% slideseq_ct_spatial_genes$gene, 1, 0)) |>
  group_by(xchr) |>
  summarise(sum(in_slideseq))
