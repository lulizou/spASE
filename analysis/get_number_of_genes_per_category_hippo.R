library(dplyr)
library(tidyr)

gencode <- import('/rafalab/lzou/resources/gencode.vM10.annotation.gff3.gz')
xchr_genes <- unique(gencode$gene_name[which(seqnames(gencode)=='chrX')])
xchr_genes <- c(xchr_genes, 'Bex3')

# Mouse 1

overall_bias_df <- readRDS('betabin_results_overall_bias_hippo_1.rds') |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))


# Total number of genes

overall_bias_df |>
  filter(totalUMI > 2^7) |>
  group_by(xchr) |>
  summarise(n())

mouse1_all_genes <- overall_bias_df |>
  filter(totalUMI > 2^7) |>
  select(gene,xchr)

# Overall p \neq 0.5

overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  group_by(xchr) |>
  summarise(n())

mouse1_bias_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  select(gene,xchr)

# Maternal bias

overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  group_by(xchr) |>
  summarise(n())

mouse1_mat_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  select(gene, xchr)

# Paternal bias

overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  group_by(xchr) |>
  summarise(n())

mouse1_pat_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  select(gene, xchr)


# Within cell type p \neq 0.5

celltype_bias_df <- readRDS('betabin_results_celltype_bias_hippo_1.rds') |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse1_ct_bias_genes <- celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Maternal cell type bias

celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse1_ct_mat_genes <- celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Paternal cell type bias

celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse1_ct_pat_genes <- celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Overall spatial pattern

overall_spatial_df <- readRDS('spasev1_results_global_combined_hippo_1.rds')$result |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

overall_spatial_df |>
  filter(qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse1_overall_spatial_genes <- overall_spatial_df |>
  filter(qval < 0.01) |>
  select(gene,xchr)

# Within cell type spatial pattern

celltype_spatial_df <- readRDS('spasev1_results_celltype_combined_hippo_1.rds')$result |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

celltype_spatial_df |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse1_ct_spatial_genes <- celltype_spatial_df |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Mouse 2

overall_bias_df <- readRDS('betabin_results_overall_bias_hippo_2.rds') |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))


# Total number of genes

overall_bias_df |>
  filter(totalUMI > 2^7) |>
  group_by(xchr) |>
  summarise(n())

mouse2_all_genes <- overall_bias_df |>
  filter(totalUMI > 2^7) |>
  select(gene,xchr)

# Overall p \neq 0.5

overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  group_by(xchr) |>
  summarise(n())

mouse2_bias_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  select(gene,xchr)

# Maternal bias

overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  group_by(xchr) |>
  summarise(n())

mouse2_mat_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  select(gene, xchr)

# Paternal bias

overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  group_by(xchr) |>
  summarise(n())

mouse2_pat_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  select(gene, xchr)


# Within cell type p \neq 0.5

celltype_bias_df <- readRDS('betabin_results_celltype_bias_hippo_2.rds') |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse2_ct_bias_genes <- celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Maternal cell type bias

celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse2_ct_mat_genes <- celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Paternal cell type bias

celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse2_ct_pat_genes <- celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)


# Overall spatial pattern

overall_spatial_df <- readRDS('spasev1_results_hippo_2.rds')$result |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

overall_spatial_df |>
  filter(qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse2_overall_spatial_genes <- overall_spatial_df |>
  filter(qval < 0.01) |>
  select(gene,xchr)

# Within cell type spatial pattern

celltype_spatial_df <- readRDS('spasev1_results_celltype_hippo_2.rds')$result |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

celltype_spatial_df |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse2_ct_spatial_genes <- celltype_spatial_df |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Mouse 3

overall_bias_df <- readRDS('betabin_results_overall_bias_hippo_3.rds') |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))


# Total number of genes

overall_bias_df |>
  filter(totalUMI > 2^7) |>
  group_by(xchr) |>
  summarise(n())

mouse3_all_genes <- overall_bias_df |>
  filter(totalUMI > 2^7) |>
  select(gene,xchr)

# Overall p \neq 0.5

overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  group_by(xchr) |>
  summarise(n())

mouse3_bias_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  select(gene,xchr)

# Maternal bias

overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  group_by(xchr) |>
  summarise(n())

mouse3_mat_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  select(gene, xchr)

# Paternal bias

overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  group_by(xchr) |>
  summarise(n())

mouse3_pat_genes <- overall_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  select(gene, xchr)


# Within cell type p \neq 0.5

celltype_bias_df <- readRDS('betabin_results_celltype_bias_hippo_3.rds') |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse3_ct_bias_genes <- celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75 | p < 0.25) ) |
           grepl('mono',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Maternal cell type bias

celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse3_ct_mat_genes <- celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p > 0.75) ) |
           grepl('monoallelic1',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)

# Paternal cell type bias

celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse3_ct_pat_genes <- celltype_bias_df |>
  filter(totalUMI > 2^7, (qval < 0.01 & (p < 0.25) ) |
           grepl('monoallelic2',flag))|>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)


# Overall spatial pattern

overall_spatial_df <- readRDS('spasev1_results_global_combined_hippo_3.rds')$result |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

overall_spatial_df |>
  filter(qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse3_overall_spatial_genes <- overall_spatial_df |>
  filter(qval < 0.01) |>
  select(gene,xchr)

# Within cell type spatial pattern

celltype_spatial_df <- readRDS('spasev1_results_celltype_combined_hippo_3.rds')$result |>
  filter(!grepl('mt-', gene)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, T, F))

celltype_spatial_df |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  group_by(xchr) |>
  summarise(n())

mouse3_ct_spatial_genes <- celltype_spatial_df |>
  group_by(gene, xchr) |>
  summarise(min_qval = min(qval)) |>
  filter(min_qval < 0.01) |>
  select(gene,xchr)


# Overlap

# Total genes

mouse1_all_genes |>
  filter(gene %in% mouse3_all_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

mouse2_all_genes |>
  filter(gene %in% mouse3_all_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

# No significant ASE

mouse1_all_genes |>
  filter(!gene %in% c(mouse1_bias_genes$gene, mouse1_ct_bias_genes$gene,
                      mouse1_overall_spatial_genes$gene, mouse1_ct_spatial_genes$gene)) |>
  group_by(xchr) |>
  summarise(n())

mouse2_all_genes |>
  filter(!gene %in% c(mouse2_bias_genes$gene, mouse2_ct_bias_genes$gene,
                      mouse2_overall_spatial_genes$gene, mouse2_ct_spatial_genes$gene)) |>
  group_by(xchr) |>
  summarise(n())

mouse3_all_genes |>
  filter(!gene %in% c(mouse3_bias_genes$gene, mouse3_ct_bias_genes$gene,
                      mouse3_overall_spatial_genes$gene, mouse3_ct_spatial_genes$gene)) |>
  group_by(xchr) |>
  summarise(n())

mouse1_no_ase_genes <- mouse1_all_genes |>
  filter(!gene %in% c(mouse1_bias_genes$gene, mouse1_ct_bias_genes$gene,
                      mouse1_overall_spatial_genes$gene, mouse1_ct_spatial_genes$gene)) |>
  select(gene, xchr)

mouse2_no_ase_genes <- mouse2_all_genes |>
  filter(!gene %in% c(mouse2_bias_genes$gene, mouse2_ct_bias_genes$gene,
                      mouse2_overall_spatial_genes$gene, mouse2_ct_spatial_genes$gene)) |>
  select(gene, xchr)

mouse3_no_ase_genes <- mouse3_all_genes |>
  filter(!gene %in% c(mouse3_bias_genes$gene, mouse3_ct_bias_genes$gene,
                      mouse3_overall_spatial_genes$gene, mouse3_ct_spatial_genes$gene)) |>
  select(gene, xchr)

mouse1_no_ase_genes |>
  filter(gene %in% mouse3_no_ase_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

mouse2_no_ase_genes |>
  filter(gene %in% mouse3_no_ase_genes$gene) |>
  group_by(xchr) |>
  summarise(n())


# Maternal bias

mouse1_mat_genes |>
  filter(gene %in% mouse3_mat_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

mouse2_mat_genes |>
  filter(gene %in% mouse3_mat_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

# Paternal bias

mouse1_pat_genes |>
  filter(gene %in% mouse3_pat_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

mouse2_pat_genes |>
  filter(gene %in% mouse3_pat_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

# Cell type bias

mouse1_ct_mat_genes |>
  filter(gene %in% mouse3_ct_mat_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

mouse2_ct_mat_genes |>
  filter(gene %in% mouse3_ct_mat_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

mouse1_ct_pat_genes |>
  filter(gene %in% mouse3_ct_pat_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

mouse2_ct_pat_genes |>
  filter(gene %in% mouse3_ct_pat_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

# Overall spatial

mouse1_overall_spatial_genes |>
  filter(gene %in% mouse3_overall_spatial_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

mouse2_overall_spatial_genes |>
  filter(gene %in% mouse3_overall_spatial_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

# Cell type spatial

mouse1_ct_spatial_genes |>
  filter(gene %in% mouse3_ct_spatial_genes$gene) |>
  group_by(xchr) |>
  summarise(n())

mouse2_ct_spatial_genes |>
  filter(gene %in% mouse3_ct_spatial_genes$gene) |>
  group_by(xchr) |>
  summarise(n())
