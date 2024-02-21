Overall results summary figures
================

- <a href="#overall-p" id="toc-overall-p">Overall p</a>
  - <a href="#estimated-p-pairs" id="toc-estimated-p-pairs">Estimated p
    pairs</a>
- <a href="#overall-p-bias" id="toc-overall-p-bias">Overall p bias</a>
- <a href="#within-cell-type-p-bias"
  id="toc-within-cell-type-p-bias">Within cell type p bias</a>
- <a href="#overall-spatial-patterns"
  id="toc-overall-spatial-patterns">Overall spatial patterns</a>
- <a href="#within-cell-type-spatial-patterns"
  id="toc-within-cell-type-spatial-patterns">Within cell type spatial
  patterns</a>
- <a href="#table-1" id="toc-table-1">Table 1</a>
- <a href="#table-2" id="toc-table-2">Table 2</a>
- <a href="#supp-table-samples" id="toc-supp-table-samples">Supp table
  samples</a>
- <a href="#supp-table-overall-bias" id="toc-supp-table-overall-bias">Supp
  table overall bias</a>
- <a href="#supp-table-ct-bias" id="toc-supp-table-ct-bias">Supp table ct
  bias</a>
- <a href="#supp-table-spatial" id="toc-supp-table-spatial">Supp table
  spatial</a>
- <a href="#supp-table-ct-spatial" id="toc-supp-table-ct-spatial">Supp
  table ct spatial</a>

``` r
library(spacexr)
library(spASE)
```

    Registered S3 method overwritten by 'spASE':
      method             from   
      merge.RCTD.objects spacexr


    Attaching package: 'spASE'

    The following objects are masked from 'package:spacexr':

        aggregate_cell_types, build.designmatrix.intercept,
        build.designmatrix.nonparam, build.designmatrix.regions,
        build.designmatrix.single, choose_sigma_c, convert.old.RCTD,
        count_cell_types, create_RCTD_plots, create.RCTD,
        create.RCTD.replicates, CSIDE.population.inference,
        exvar.celltocell.interactions, exvar.point.density, fitBulk,
        fitPixels, get_cell_type_info, get_de_genes, get_doublet_weights,
        get_norm_ref, get_standard_errors, import_weights,
        make_all_de_plots, make_de_plots_genes, make_de_plots_quant,
        make_de_plots_regions, make_de_plots_replicates,
        make_de_plots_spatial, normalize_weights, plot_all_cell_types,
        plot_class, plot_cond_occur, plot_doub_occur_stack, plot_doublets,
        plot_doublets_type, plot_gene_raw, plot_gene_regions,
        plot_gene_two_regions, plot_occur_unthreshold,
        plot_prediction_gene, plot_puck_continuous, plot_puck_wrapper,
        plot_weights, plot_weights_doublet, plot_weights_unthreshold,
        process_beads_batch, process_data, read.SpatialRNA,
        read.VisiumSpatialRNA, Reference, restrict_counts, restrict_puck,
        run.CSIDE, run.CSIDE.general, run.CSIDE.intercept,
        run.CSIDE.nonparam, run.CSIDE.regions, run.CSIDE.replicates,
        run.CSIDE.single, run.RCTD, run.RCTD.replicates,
        save.CSIDE.replicates, set_cell_types_assigned,
        set_likelihood_vars, SpatialRNA, write_de_summary

``` r
library(dplyr)
```


    Attaching package: 'dplyr'

    The following objects are masked from 'package:stats':

        filter, lag

    The following objects are masked from 'package:base':

        intersect, setdiff, setequal, union

``` r
library(ggplot2)
library(data.table)
```


    Attaching package: 'data.table'

    The following objects are masked from 'package:dplyr':

        between, first, last

``` r
library(viridis)
```

    Loading required package: viridisLite

``` r
library(gplots)
```


    Attaching package: 'gplots'

    The following object is masked from 'package:stats':

        lowess

``` r
library(rtracklayer)
```

    Loading required package: GenomicRanges

    Loading required package: stats4

    Loading required package: BiocGenerics


    Attaching package: 'BiocGenerics'

    The following objects are masked from 'package:dplyr':

        combine, intersect, setdiff, union

    The following objects are masked from 'package:stats':

        IQR, mad, sd, var, xtabs

    The following objects are masked from 'package:base':

        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
        table, tapply, union, unique, unsplit, which.max, which.min

    Loading required package: S4Vectors


    Attaching package: 'S4Vectors'

    The following object is masked from 'package:gplots':

        space

    The following objects are masked from 'package:data.table':

        first, second

    The following objects are masked from 'package:dplyr':

        first, rename

    The following object is masked from 'package:utils':

        findMatches

    The following objects are masked from 'package:base':

        expand.grid, I, unname

    Loading required package: IRanges


    Attaching package: 'IRanges'

    The following object is masked from 'package:data.table':

        shift

    The following objects are masked from 'package:dplyr':

        collapse, desc, slice

    Loading required package: GenomeInfoDb

``` r
library(Matrix)
```


    Attaching package: 'Matrix'

    The following object is masked from 'package:S4Vectors':

        expand

``` r
THRESHOLD = 500
```

``` r
# Read in gencode to grab xchr genes
gencode <- import('results/gencode.vM10.annotation.gff3.gz')
xchr_genes <- unique(gencode$gene_name[which(seqnames(gencode)=='chrX')])
xchr_genes <- c(xchr_genes, 'Bex3')
```

# Overall p

``` r
hippo_1 <- readRDS('results/results_overall_bias_hippo_1.rds')
hippo_2 <- readRDS('results/results_overall_bias_hippo_2.rds')
hippo_3 <- readRDS('results/results_overall_bias_hippo_3.rds')
cere_3 <- readRDS('results/results_overall_bias_cere_3.rds')
cere_4 <- readRDS('results/results_overall_bias_cere_4_visium.rds')
```

``` r
df <- hippo_1 |> 
  filter(totalUMI > THRESHOLD) |>
  mutate(hippo_1 = p) |>
  select(gene, hippo_1) |>
  left_join(hippo_2 |>
              filter(totalUMI > THRESHOLD) |>
              mutate(hippo_2 = p) |>
              select(gene, hippo_2)) |>
  left_join(hippo_3 |>
              filter(totalUMI > THRESHOLD) |>
              mutate(hippo_3 = p) |>
              select(gene, hippo_3)) |>
  left_join(cere_3 |>
              filter(totalUMI > THRESHOLD) |>
              mutate(cere_3 = p) |>
              select(gene, cere_3)) |>
  left_join(cere_4 |>
              filter(totalUMI > THRESHOLD) |>
              mutate(cere_4 = p) |>
              select(gene, cere_4)) |>
  filter(!is.na(hippo_1), !is.na(hippo_2), !is.na(hippo_3), !is.na(cere_3), !is.na(cere_4)) |>
  mutate(xchr = ifelse(gene %in% xchr_genes, 1, 0)) |>
  filter(!grepl('mt-', gene))
```

    Joining with `by = join_by(gene)`
    Joining with `by = join_by(gene)`
    Joining with `by = join_by(gene)`
    Joining with `by = join_by(gene)`

``` r
mm <- as.matrix(df[,2:6])
rownames(mm) <- df$gene
colnames(mm) <- c('Hippo 1', 'Hippo 2', 'Hippo 3', 'Cere 3', 'Cere 4')
pdf('figures/04_heatmap.pdf', height=8, width=4)
heatmap.2(mm, scale='none', col= bluered(100), trace='none', density.info='none', breaks = seq(0,1,length.out=101), lhei=c(1,6), cexRow = 0.5, cexCol = 0.5)
dev.off()
```

    quartz_off_screen 
                    2 

## Estimated p pairs

``` r
panel.hist <- function(x, ...)
{
  usr <- par("usr")
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE, breaks = 50)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))^2
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor)
}
pdf('figures/04_pair_scatter.pdf', height=6, width=6)
pairs(df |> filter(xchr == 0) |> select(-gene, -xchr) |> dplyr::rename(`Hippo 1` = hippo_1, `Hippo 2` = hippo_2, `Hippo 3` = hippo_3, `Cere 3` = cere_3, `Cere 4` = cere_4) |> as.matrix(),
      upper.panel = panel.cor,
      diag.panel = panel.hist)
dev.off()
```

    quartz_off_screen 
                    2 

# Overall p bias

``` r
extract_mat_pat_overall <- function(spase_res, qthresh = 0.01) {
  dd <- spase_res |> 
    filter(qval < qthresh | grepl('monoallelic', flag)) |>
    mutate(bias = case_when(
      p > 0.6 ~ 'maternal',
      p < 0.4 ~ 'paternal',
      flag == 'monoallelic1' ~ 'maternal',
      flag == 'monoallelic2' ~ 'paternal',
      T ~ 'low_bias'
    ))
  dd$xchr <- ifelse(dd$gene %in% xchr_genes, T, F)
  return(dd)
}
hippo1_overall <- extract_mat_pat_overall(hippo_1)

hippo2_overall <- extract_mat_pat_overall(hippo_2)

hippo3_overall <- extract_mat_pat_overall(hippo_3)

cere4_overall <- extract_mat_pat_overall(cere_4)

cere3_overall <- extract_mat_pat_overall(cere_3)
```

# Within cell type p bias

``` r
hippo1_intercept <- readRDS('results/results_celltype_intercept_hippo_1_df_5.rds')
hippo2_intercept <- readRDS('results/results_celltype_intercept_hippo_2_df_5.rds')
hippo3_intercept <- readRDS('results/results_celltype_intercept_hippo_3_df_5.rds')
cere3_intercept <- readRDS('results/results_celltype_intercept_cere_3_df_5.rds')
cere4_intercept <- readRDS('results/results_celltype_intercept_cere_4_visium_df_5.rds')
```

``` r
extract_ct_genes <- function(sig_gene_list) {
  return(unique(unname(unlist(lapply(sig_gene_list, rownames)))))
}

extract_mat_pat_ct <- function(sig_gene_list) {
  res <- NULL
  for (i in seq_along(sig_gene_list)) {
    dd <- sig_gene_list[[i]]
    if (nrow(dd)==0) {
      next
    }
    ct <- names(sig_gene_list)[i]
    dd$gene <- rownames(dd)
    rownames(dd) <- NULL
    dd$cell_type <- ct
    dd$bias <- ifelse(dd$log_fc > 0, 'maternal', 'paternal')
    dd$xchr <- ifelse(dd$gene %in% xchr_genes, T, F)
    if (is.null(res)) {
      res <- dd
    } else {
      res <- bind_rows(res, dd)
    }
  }
  res <- res |> select(gene, xchr, cell_type, bias, Z_score, log_fc, se, p_val, mean_0, mean_1, sd_0, sd_1)
  return(res)
}
```

``` r
hippo1_ct <- extract_mat_pat_ct(hippo1_intercept@spase_results$sig_gene_list)

hippo2_ct <- extract_mat_pat_ct(hippo2_intercept@spase_results$sig_gene_list)

hippo3_ct <- extract_mat_pat_ct(hippo3_intercept@spase_results$sig_gene_list)

cere3_ct <- extract_mat_pat_ct(cere3_intercept@spase_results$sig_gene_list)

cere4_ct <- extract_mat_pat_ct(cere4_intercept@spase_results$sig_gene_list)

rm(hippo1_intercept, hippo2_intercept, hippo3_intercept, cere3_intercept, cere4_intercept)
```

# Overall spatial patterns

``` r
hippo_1_sp <- readRDS('results/results_overall_spatial_hippo_1.rds')
hippo_2_sp <- readRDS('results/results_overall_spatial_hippo_2.rds')
hippo_3_sp <- readRDS('results/results_overall_spatial_hippo_3.rds')
cere_3_sp <- readRDS('results/results_overall_spatial_cere_3.rds')
cere_4_sp <- readRDS('results/results_overall_spatial_cere_4_visium.rds')
```

``` r
extract_spatial_overall <- function(spase_res, qthresh = 0.01) {
  dd <- spase_res |> filter(qval < qthresh)
  dd$xchr <- ifelse(dd$gene %in% xchr_genes, T, F)
  return(dd)
}

hippo1_sp_overall <- extract_spatial_overall(hippo_1_sp$result)

hippo2_sp_overall <- extract_spatial_overall(hippo_2_sp$result)

hippo3_sp_overall <- extract_spatial_overall(hippo_3_sp$result)

cere3_sp_overall <- extract_spatial_overall(cere_3_sp$result)

cere4_sp_overall <- extract_spatial_overall(cere_4_sp$result)
```

# Within cell type spatial patterns

``` r
hippo_1_spct_res <- readRDS('results/results_celltype_hippo_1_df_5.rds')
hippo_2_spct_res <- readRDS('results/results_celltype_hippo_2_df_5.rds')
hippo_3_spct_res <- readRDS('results/results_celltype_hippo_3_df_5.rds')
cere_3_spct_res <- readRDS('results/results_celltype_cere_3_df_5.rds')
cere_4_spct_res <- readRDS('results/results_celltype_cere_4_visium_df_5.rds')
```

``` r
extract_spatial_ct <- function(sig_gene_list) {
  res <- NULL
  for (i in seq_along(sig_gene_list)) {
    dd <- sig_gene_list[[i]]
    if (nrow(dd) == 0) {
      next
    }
    ct <- names(sig_gene_list)[i]
    dd$gene <- rownames(dd)
    rownames(dd) <- NULL
    dd$cell_type <- ct
    dd$xchr <- ifelse(dd$gene %in% xchr_genes, T, F)
    if (is.null(res)) {
      res <- dd
    } else {
      res <- bind_rows(res, dd)
    }
  }
  return(res)
}


hippo1_sp_ct <- extract_spatial_ct(hippo_1_spct_res@spase_results$sig_gene_list)

hippo2_sp_ct <- extract_spatial_ct(hippo_2_spct_res@spase_results$sig_gene_list)

hippo3_sp_ct <- extract_spatial_ct(hippo_3_spct_res@spase_results$sig_gene_list)

cere3_sp_ct <- extract_spatial_ct(cere_3_spct_res@spase_results$sig_gene_list)

cere4_sp_ct <- extract_spatial_ct(cere_4_spct_res@spase_results$sig_gene_list)

hippo1_sp_ct_all <- extract_spatial_ct(hippo_1_spct_res@spase_results$all_gene_list)

hippo2_sp_ct_all <- extract_spatial_ct(hippo_2_spct_res@spase_results$all_gene_list)

hippo3_sp_ct_all <- extract_spatial_ct(hippo_3_spct_res@spase_results$all_gene_list)

cere3_sp_ct_all <- extract_spatial_ct(cere_3_spct_res@spase_results$all_gene_list)

cere4_sp_ct_all <- extract_spatial_ct(cere_4_spct_res@spase_results$all_gene_list)
```

``` r
format_lengths <- function(c1, c2) {
  return(paste0(prettyNum(length(c1), big.mark=','), ' (', prettyNum(length(c2), big.mark=','), ')'))
}
extract_info_for_table <- function(all, p_bias, ct_p_bias, sp, sp_ct, sp_ct_all) {
  all$xchr <- ifelse(all$gene %in% xchr_genes, T, F)
  all_genes_autosome <- all |> filter(!xchr) |> pull(gene)
  rn <- rownames(sp_ct_all@originalSpatialRNA@maternalCounts[rowSums(sp_ct_all@originalSpatialRNA@maternalCounts) + rowSums(sp_ct_all@originalSpatialRNA@paternalCounts) > 128,])
  rn_autosome <- rn[!rn %in% xchr_genes]
  rn_xchr <- rn[rn %in% xchr_genes]
  all_genes_autosome <- unique(c(all_genes_autosome, rn_autosome))
  all_genes_xchr <- all |> filter(xchr) |> pull(gene)
  all_genes_xchr <- unique(c(all_genes_xchr, rn_xchr))
  mat_genes_autosome <- p_bias |> filter(!xchr, bias == 'maternal') |> pull(gene)
  pat_genes_autosome <- p_bias |> filter(!xchr, bias == 'paternal') |> pull(gene)
  mat_genes_xchr <- p_bias |> filter(xchr, bias == 'maternal') |> pull(gene)
  pat_genes_xchr <- p_bias |> filter(xchr, bias == 'paternal') |> pull(gene)
  mat_genes_ct_autosome <- ct_p_bias |> filter(!xchr, bias == 'maternal') |> pull(gene) |> unique()
  pat_genes_ct_autosome <- ct_p_bias |> filter(!xchr, bias == 'paternal') |> pull(gene)|> unique()
  mat_genes_ct_xchr <- ct_p_bias |> filter(xchr, bias == 'maternal') |> pull(gene)|> unique()
  pat_genes_ct_xchr <- ct_p_bias |> filter(xchr, bias == 'paternal') |> pull(gene)|> unique()
  all_spatial_autosome <- sp |> filter(!xchr) |> pull(gene)
  all_spatial_xchr <- sp |> filter(xchr) |> pull(gene)
  if (is.null(sp_ct)) {
    ct_spatial_autosome <- c()
    ct_spatial_xchr <- c()
  } else {
    ct_spatial_autosome <- sp_ct |> filter(!xchr) |> pull(gene) |> unique()
    ct_spatial_xchr <- sp_ct |> filter(xchr) |> pull(gene) |> unique()
  }
  no_sig_ase_autosome <- all_genes_autosome[!all_genes_autosome %in% c(mat_genes_autosome, pat_genes_autosome, mat_genes_ct_autosome, pat_genes_ct_autosome, all_spatial_autosome, ct_spatial_autosome)]
  no_sig_ase_xchr <- all_genes_xchr[!all_genes_xchr %in% c(mat_genes_xchr, pat_genes_xchr, mat_genes_ct_xchr, pat_genes_ct_xchr, all_spatial_xchr, ct_spatial_xchr)]
  return(list(
    col = c(
      format_lengths(no_sig_ase_autosome, no_sig_ase_xchr),
      format_lengths(mat_genes_autosome, mat_genes_xchr),
      format_lengths(pat_genes_autosome, pat_genes_xchr),
      format_lengths(mat_genes_ct_autosome, mat_genes_ct_xchr),
      format_lengths(pat_genes_ct_autosome, pat_genes_ct_xchr),
      format_lengths(all_spatial_autosome, all_spatial_xchr),
      format_lengths(ct_spatial_autosome, ct_spatial_xchr),
      format_lengths(all_genes_autosome, all_genes_xchr)
    ),
    nsa = no_sig_ase_autosome, nsx = no_sig_ase_xchr,
    mga = mat_genes_autosome, mgx = mat_genes_xchr, pga = pat_genes_autosome, pgx = pat_genes_xchr, mgca = mat_genes_ct_autosome, mgcx = mat_genes_ct_xchr, pgca = pat_genes_ct_autosome, pgcx = pat_genes_ct_xchr, asa = all_spatial_autosome, asx = all_spatial_xchr, cta = ct_spatial_autosome, ctx = ct_spatial_xchr,
    aga = all_genes_autosome, agx = all_genes_xchr
  ))
}
```

``` r
h1 <- extract_info_for_table(hippo_1, hippo1_overall, hippo1_ct, hippo1_sp_overall, hippo1_sp_ct, hippo_1_spct_res)
h2 <- extract_info_for_table(hippo_2, hippo2_overall, hippo2_ct, hippo2_sp_overall, hippo2_sp_ct, hippo_2_spct_res)
h3 <- extract_info_for_table(hippo_3, hippo3_overall, hippo3_ct, hippo3_sp_overall, hippo3_sp_ct, hippo_3_spct_res)
c3 <- extract_info_for_table(cere_3, cere3_overall, cere3_ct, cere3_sp_overall, cere3_sp_ct, cere_3_spct_res)
c4 <- extract_info_for_table(cere_4, cere4_overall, cere4_ct, cere4_sp_overall, cere4_sp_ct, cere_4_spct_res)

# make the overlap columns
c_overlap <- c(
  format_lengths(c3$nsa[c3$nsa %in% c4$nsa], c3$nsx[c3$nsx %in% c4$nsx]),
  format_lengths(c3$mga[c3$mga %in% c4$mga], c3$mgx[c3$mgx %in% c4$mgx]),
  format_lengths(c3$pga[c3$pga %in% c4$pga], c3$pgx[c3$pgx %in% c4$pgx]),
  format_lengths(c3$mgca[c3$mgca %in% c4$mgca], c3$mgcx[c3$mgcx %in% c4$mgcx]),
  format_lengths(c3$pgca[c3$pgca %in% c4$pgca], c3$pgcx[c3$pgcx %in% c4$pgcx]),
  format_lengths(c3$asa[c3$asa %in% c4$asa], c3$asx[c3$asx %in% c4$asx]),
  format_lengths(c3$cta[c3$cta %in% c4$cta], c3$ctx[c3$ctx %in% c4$ctx]),
  format_lengths(c3$aga[c3$aga %in% c4$aga], c3$agx[c3$agx %in% c4$agx])
)

h_overlap12 <- c(
  format_lengths(h1$nsa[h1$nsa %in% h2$nsa], h1$nsx[h1$nsx %in% h2$nsx]),
  format_lengths(h1$mga[h1$mga %in% h2$mga], h1$mgx[h1$mgx %in% h2$mgx]),
  format_lengths(h1$pga[h1$pga %in% h2$pga], h1$pgx[h1$pgx %in% h2$pgx]),
  format_lengths(h1$mgca[h1$mgca %in% h2$mgca], h1$mgcx[h1$mgcx %in% h2$mgcx]),
  format_lengths(h1$pgca[h1$pgca %in% h2$pgca], h1$pgcx[h1$pgcx %in% h2$pgcx]),
  format_lengths(h1$asa[h1$asa %in% h2$asa], h1$asx[h1$asx %in% h2$asx]),
  format_lengths(h1$cta[h1$cta %in% h2$cta], h1$ctx[h1$ctx %in% h2$ctx]),
  format_lengths(h1$aga[h1$aga %in% h2$aga], h1$agx[h1$agx %in% h2$agx])
)

h_overlap13 <- c(
  format_lengths(h1$nsa[h1$nsa %in% h3$nsa], h1$nsx[h1$nsx %in% h3$nsx]),
  format_lengths(h1$mga[h1$mga %in% h3$mga], h1$mgx[h1$mgx %in% h3$mgx]),
  format_lengths(h1$pga[h1$pga %in% h3$pga], h1$pgx[h1$pgx %in% h3$pgx]),
  format_lengths(h1$mgca[h1$mgca %in% h3$mgca], h1$mgcx[ h1$mgcx %in% h3$mgcx]),
  format_lengths(h1$pgca[h1$pgca %in% h3$pgca], h1$pgcx[h1$pgcx %in% h3$pgcx]),
  format_lengths(h1$asa[h1$asa %in% h3$asa], h1$asx[h1$asx %in% h3$asx]),
  format_lengths(h1$cta[h1$cta %in% h3$cta], h1$ctx[h1$ctx %in% h3$ctx]),
  format_lengths(h1$aga[h1$aga %in% h3$aga], h1$agx[h1$agx %in% h3$agx])
)
```

# Table 1

``` r
t1 <- data.frame(Category = c('No significant ASE', 'Overall maternal bias', 'Overall paternal bias', 'Within cell type maternal bias', 'Within cell type paternal bias', 'Overall spatial pattern', 'Within cell type spatial pattern', 'Total n genes'),
                 `Slide-seq (Mouse 3)` = c3$col, `Visium (Mouse 4)` = c4$col, Overlap = c_overlap)
print(xtable::xtable(t1))
```

    % latex table generated in R 4.3.1 by xtable 1.8-4 package
    % Wed Feb 21 11:43:57 2024
    \begin{table}[ht]
    \centering
    \begin{tabular}{rllll}
      \hline
     & Category & Slide.seq..Mouse.3. & Visium..Mouse.4. & Overlap \\ 
      \hline
    1 & No significant ASE & 6,371 (55) & 8,972 (162) & 5,544 (31) \\ 
      2 & Overall maternal bias & 720 (157) & 502 (112) & 196 (81) \\ 
      3 & Overall paternal bias & 947 (7) & 672 (12) & 300 (5) \\ 
      4 & Within cell type maternal bias & 216 (43) & 4 (1) & 3 (1) \\ 
      5 & Within cell type paternal bias & 290 (3) & 8 (0) & 3 (0) \\ 
      6 & Overall spatial pattern & 8 (19) & 2 (0) & 1 (0) \\ 
      7 & Within cell type spatial pattern & 0 (5) & 0 (0) & 0 (0) \\ 
      8 & Total n genes & 8,304 (225) & 10,147 (286) & 7,599 (204) \\ 
       \hline
    \end{tabular}
    \end{table}

# Table 2

``` r
t2 <- data.frame(Category = c('No significant ASE', 'Overall maternal bias', 'Overall paternal bias', 'Within cell type maternal bias', 'Within cell type paternal bias', 'Overall spatial pattern', 'Within cell type spatial pattern', 'Total n genes'),
                 `Mouse 1` = h1$col, `Mouse 2` = h2$col, `Mouse 3` = h3$col, Overlap12 = h_overlap12, Overlap13 = h_overlap13)
print(xtable::xtable(t2))
```

    % latex table generated in R 4.3.1 by xtable 1.8-4 package
    % Wed Feb 21 11:43:57 2024
    \begin{table}[ht]
    \centering
    \begin{tabular}{rllllll}
      \hline
     & Category & Mouse.1 & Mouse.2 & Mouse.3 & Overlap12 & Overlap13 \\ 
      \hline
    1 & No significant ASE & 4,965 (104) & 2,242 (5) & 6,549 (16) & 2,012 (3) & 4,264 (4) \\ 
      2 & Overall maternal bias & 349 (30) & 176 (55) & 623 (208) & 93 (11) & 179 (29) \\ 
      3 & Overall paternal bias & 456 (9) & 180 (1) & 834 (1) & 109 (1) & 254 (1) \\ 
      4 & Within cell type maternal bias & 101 (2) & 23 (14) & 229 (63) & 17 (0) & 61 (1) \\ 
      5 & Within cell type paternal bias & 146 (9) & 21 (0) & 307 (0) & 18 (0) & 94 (0) \\ 
      6 & Overall spatial pattern & 18 (17) & 1 (0) & 67 (0) & 1 (0) & 9 (0) \\ 
      7 & Within cell type spatial pattern & 0 (0) & 0 (0) & 1 (0) & 0 (0) & 0 (0) \\ 
      8 & Total n genes & 5,866 (159) & 2,609 (61) & 8,309 (225) & 2,560 (59) & 5,698 (150) \\ 
       \hline
    \end{tabular}
    \end{table}

Note that is detected as significant within Interneuron for Mouse 3;
however, this goes away when controlling for the Interneuron sub-type
located in the cluster of cells with high Sst expression. Thus, it is
not included in the final table.

# Supp table samples

``` r
get_stats <- function(cside_spase_obj) {
  tot <- cside_spase_obj@originalSpatialRNA@counts
  tot_ase <- cside_spase_obj@originalSpatialRNA@maternalCounts + cside_spase_obj@originalSpatialRNA@paternalCounts
  nspots <- ncol(tot)
  readsspot <- mean(colSums(tot))
  readsase <- mean(colSums(tot_ase))
  return(c(prettyNum(nspots, big.mark=','), prettyNum(round(readsspot),big.mark=','), prettyNum(round(readsase),big.mark=',')))
}

rr <- rbind(get_stats(hippo_1_spct_res), get_stats(hippo_2_spct_res), get_stats(hippo_3_spct_res), get_stats(cere_3_spct_res), get_stats(cere_4_spct_res))
colnames(rr) <- c('N spots', 'Avg. reads/spot', 'Allele-resolved')
rr <- data.frame(rr)

dd <- bind_cols(
  data.frame(Mouse = c(1,2,3,3,4), 
             Tissue = c('Hippocampus', 'Hippocampus', 'Hippocampus', 'Cerebellum', 'Cerebellum'), 
             Platform = c('Slide-seq', 'Slide-seq', 'Slide-seq', 'Slide-seq', 'Visium'), 
             `Read length` = c(160, 56, 250, 250, 91)),
  rr
)
xtable::xtable(dd)
```

    % latex table generated in R 4.3.1 by xtable 1.8-4 package
    % Wed Feb 21 11:44:04 2024
    \begin{table}[ht]
    \centering
    \begin{tabular}{rrllrlll}
      \hline
     & Mouse & Tissue & Platform & Read.length & N.spots & Avg..reads.spot & Allele.resolved \\ 
      \hline
    1 & 1.00 & Hippocampus & Slide-seq & 160.00 & 26,429 & 459 & 207 \\ 
      2 & 2.00 & Hippocampus & Slide-seq & 56.00 & 13,680 & 623 & 145 \\ 
      3 & 3.00 & Hippocampus & Slide-seq & 250.00 & 78,806 & 552 & 218 \\ 
      4 & 3.00 & Cerebellum & Slide-seq & 250.00 & 60,942 & 622 & 240 \\ 
      5 & 4.00 & Cerebellum & Visium & 91.00 & 4,315 & 2,219 & 834 \\ 
       \hline
    \end{tabular}
    \end{table}

# Supp table overall bias

``` r
st1 <- hippo1_overall |>
  mutate(tissue = 'Hippocampus', sample = 'Mouse 1', technology = 'Slide-seq') |>
  select(tissue, sample, technology, gene, totalUMI, totalCells, p, ci.low, ci.high, z, pval, qval, bias, xchr) |>
  bind_rows(
    hippo2_overall |>
      mutate(tissue = 'Hippocampus', sample = 'Mouse 2', technology = 'Slide-seq') |>
      select(tissue, sample, technology, gene, totalUMI, totalCells, p, ci.low, ci.high, z, pval, qval, bias, xchr)
  ) |>
  bind_rows(
    hippo3_overall |>
      mutate(tissue = 'Hippocampus', sample = 'Mouse 3', technology = 'Slide-seq') |>
      select(tissue, sample, technology, gene, totalUMI, totalCells, p, ci.low, ci.high, z, pval, qval, bias, xchr)
  ) |>
  bind_rows(
    cere3_overall |>
      mutate(tissue = 'Cerebellum', sample = 'Mouse 3', technology = 'Slide-seq') |>
      select(tissue, sample, technology, gene, totalUMI, totalCells, p, ci.low, ci.high, z, pval, qval, bias, xchr)
  ) |>
  bind_rows(
    cere4_overall |>
      mutate(tissue = 'Cerebellum', sample = 'Mouse 4', technology = 'Visium') |>
      select(tissue, sample, technology, gene, totalUMI, totalCells, p, ci.low, ci.high, z, pval, qval, bias, xchr)
  )
write.csv(st1, file = 'tables/04_supp_table_overall_bias.csv')
```

# Supp table ct bias

``` r
st2 <- hippo1_ct |>
  mutate(tissue = 'Hippocampus', sample = 'Mouse 1', technology = 'Slide-seq') |>
  select(tissue, sample, technology, gene, cell_type, Z_score, log_fc, se, p_val, bias, xchr) |>
  bind_rows(
    hippo2_ct |>
      mutate(tissue = 'Hippocampus', sample = 'Mouse 2', technology = 'Slide-seq') |>
      select(tissue, sample, technology, gene, cell_type, Z_score, log_fc, se, p_val, bias, xchr)
  ) |>
  bind_rows(
    hippo3_ct |>
      mutate(tissue = 'Hippocampus', sample = 'Mouse 3', technology = 'Slide-seq') |>
      select(tissue, sample, technology, gene, cell_type, Z_score, log_fc, se, p_val, bias, xchr
      )
  ) |>
  bind_rows(
    cere3_ct |>
      mutate(tissue = 'Cerebellum', sample = 'Mouse 3', technology = 'Slide-seq') |>
      select(tissue, sample, technology, gene, cell_type, Z_score, log_fc, se, p_val, bias, xchr)
  )|>
  bind_rows(
    cere4_ct |>
      mutate(tissue = 'Cerebellum', sample = 'Mouse 4', technology = 'Visium') |>
      select(tissue, sample, technology, gene, cell_type, Z_score, log_fc, se, p_val, bias, xchr)
  )
write.csv(st2, file = 'tables/04_supp_table_celltype_bias.csv')
```

# Supp table spatial

``` r
st3 <- hippo1_sp_overall |>
  mutate(tissue = 'Hippocampus', sample = 'Mouse 1', technology = 'Slide-seq') |>
  select(tissue, sample, technology, gene, totalUMI, totalSpots, chisq.p, qval,  xchr) |>
  bind_rows(
    hippo2_sp_overall |>
      mutate(tissue = 'Hippocampus', sample = 'Mouse 2', technology = 'Slide-seq') |>
      select(tissue, sample, technology, gene, totalUMI, totalSpots, chisq.p, qval,  xchr)
  ) |>
  bind_rows(
    hippo3_sp_overall |>
      mutate(tissue = 'Hippocampus', sample = 'Mouse 3', technology = 'Slide-seq') |>
      select(tissue, sample, technology, gene, totalUMI, totalSpots, chisq.p, qval,  xchr)
  ) |>
  bind_rows(
    cere3_sp_overall |>
      mutate(tissue = 'Cerebellum', sample = 'Mouse 3', technology = 'Slide-seq') |>
      select(tissue, sample, technology, gene, totalUMI, totalSpots, chisq.p, qval,  xchr)
  ) |>
  bind_rows(
    cere4_sp_overall |>
      mutate(tissue = 'Cerebellum', sample = 'Mouse 4', technology = 'Visium') |>
      select(tissue, sample, technology, gene, totalUMI, totalSpots, chisq.p, qval,  xchr)
  )
write.csv(st3, file = 'tables/04_supp_table_spatial_overall_bias.csv')
```

``` r
print(xtable::xtable(st3 |> select(-chisq.p) |> arrange(sample, tissue, qval), display = c('s', 's', 's', 's', 's', 'd', 'd', 'E', 's')), include.rownames=F)
```

    % latex table generated in R 4.3.1 by xtable 1.8-4 package
    % Wed Feb 21 11:44:05 2024
    \begin{table}[ht]
    \centering
    \begin{tabular}{llllrrrl}
      \hline
    tissue & sample & technology & gene & totalUMI & totalSpots & qval & xchr \\ 
      \hline
    Hippocampus & Mouse 1 & Slide-seq & Tspan7 & 9612 & 6667 & 0.00E+00 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Nrip3 & 4213 & 3060 & 9.29E-12 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Atrx & 2351 & 1893 & 1.67E-11 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Sst & 1429 & 745 & 1.67E-11 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Ptgds & 2374 & 1321 & 4.01E-10 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Arhgef9 & 1857 & 1556 & 1.09E-09 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Ids & 1002 & 919 & 1.09E-09 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Rgs4 & 2712 & 2097 & 1.54E-09 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Chgb & 7898 & 3923 & 5.53E-08 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Gprasp1 & 2665 & 2177 & 6.03E-08 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Cpe & 22568 & 11861 & 7.72E-08 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Nrsn1 & 5220 & 4092 & 1.60E-06 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Pcsk1n & 3208 & 2851 & 2.27E-06 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Maged1 & 1364 & 1235 & 3.04E-06 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Fam193a & 558 & 533 & 6.47E-06 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Crym & 1241 & 940 & 4.61E-05 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Gpm6b & 2203 & 1963 & 9.37E-05 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Sparcl1 & 18578 & 10638 & 9.37E-05 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Meg3 & 18494 & 5727 & 1.47E-04 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Bex2 & 1570 & 1391 & 2.41E-04 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Magee1 & 1698 & 1444 & 2.61E-04 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Gap43 & 2920 & 2282 & 4.71E-04 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Map1b & 14912 & 7333 & 5.28E-04 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Armcx3 & 646 & 600 & 7.98E-04 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Malat1 & 72124 & 12172 & 1.49E-03 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Rprm & 1107 & 796 & 1.49E-03 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Prps1 & 879 & 790 & 1.50E-03 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Tceal5 & 2117 & 1781 & 1.50E-03 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Frrs1l & 2992 & 2321 & 5.80E-03 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Mageh1 & 673 & 634 & 5.80E-03 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Xist & 537 & 449 & 6.37E-03 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Atp1a2 & 8717 & 5701 & 6.76E-03 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Fgf13 & 954 & 869 & 6.76E-03 & TRUE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Napa & 1299 & 1184 & 6.76E-03 & FALSE \\ 
      Hippocampus & Mouse 1 & Slide-seq & Tceal3 & 1491 & 1271 & 7.25E-03 & TRUE \\ 
      Hippocampus & Mouse 2 & Slide-seq & Ptgds & 910 & 635 & 2.96E-05 & FALSE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Meg3 & 28544 & 14574 & 0.00E+00 & FALSE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Ptgds & 6457 & 3277 & 0.00E+00 & FALSE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Rps4x & 12261 & 9019 & 3.32E-13 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Bex2 & 9202 & 7044 & 2.70E-09 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Tceal3 & 9503 & 6693 & 2.03E-08 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Cdr1 & 2255 & 1945 & 4.74E-08 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Tspan7 & 11409 & 8731 & 4.74E-08 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Bex3 & 12213 & 8369 & 7.34E-08 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Ncald & 2311 & 1804 & 1.92E-07 & FALSE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Htatsf1 & 1277 & 1068 & 1.64E-05 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Ndufb11 & 13950 & 9732 & 9.11E-05 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Tceal5 & 5929 & 4637 & 9.11E-05 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Uba1 & 7773 & 5728 & 1.73E-04 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Nkap & 962 & 783 & 2.35E-04 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Idh3g & 7591 & 5758 & 3.97E-04 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Gprasp1 & 3373 & 2746 & 7.26E-04 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Atrx & 4767 & 3617 & 7.33E-04 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Lpar1 & 1311 & 1021 & 9.02E-04 & FALSE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Mageh1 & 3337 & 2566 & 1.13E-03 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Xist & 1678 & 1373 & 1.22E-03 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Cdk16 & 3542 & 3133 & 1.22E-03 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Psat1 & 1507 & 1237 & 1.82E-03 & FALSE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Aldoc & 89897 & 26422 & 2.16E-03 & FALSE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Plp1 & 17084 & 5957 & 3.61E-03 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Shank1 & 3427 & 2879 & 3.61E-03 & FALSE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Hprt & 2848 & 2429 & 7.30E-03 & TRUE \\ 
      Cerebellum & Mouse 3 & Slide-seq & Kcnj9 & 3493 & 2655 & 7.30E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Chgb & 17790 & 8654 & 0.00E+00 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Nrip3 & 12045 & 8028 & 0.00E+00 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Ptgds & 4863 & 3318 & 0.00E+00 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Rasgrf1 & 4462 & 3399 & 0.00E+00 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Snca & 20846 & 10276 & 0.00E+00 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Sst & 6781 & 3676 & 0.00E+00 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Stmn1 & 68985 & 30606 & 0.00E+00 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & A830018L16Rik & 1831 & 1483 & 9.99E-13 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Uchl1 & 37040 & 20455 & 2.05E-12 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Lypd1 & 2073 & 1585 & 3.00E-12 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Pcp4l1 & 13914 & 9034 & 7.69E-11 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Akap5 & 196564 & 62287 & 3.71E-10 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Nrsn1 & 10799 & 8182 & 3.88E-10 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Gstm7 & 1150 & 1018 & 7.29E-10 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Mrps11 & 2140 & 1715 & 3.25E-09 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Calb1 & 1604 & 1327 & 2.49E-08 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Grb14 & 2591 & 2119 & 7.86E-08 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Sparcl1 & 40302 & 25100 & 1.21E-07 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Spon1 & 1941 & 1554 & 1.36E-07 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Bok & 2728 & 2336 & 1.03E-06 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Car12 & 502 & 396 & 2.36E-06 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Rasgrp1 & 11580 & 7368 & 2.36E-06 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Selenop & 21213 & 14535 & 2.36E-06 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Rhpn2 & 862 & 782 & 2.76E-06 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Ccng1 & 2158 & 1806 & 5.14E-06 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Gaa & 4307 & 3332 & 1.01E-05 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Cck & 15170 & 10001 & 1.34E-05 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Chst2 & 8951 & 6587 & 1.43E-05 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Rprm & 2537 & 1678 & 1.50E-05 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Nefl & 17453 & 10783 & 1.66E-05 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Actr2 & 8748 & 6525 & 2.05E-05 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Cnih3 & 733 & 641 & 3.18E-05 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Ncald & 16168 & 9706 & 6.65E-05 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Cpe & 54353 & 30714 & 7.63E-05 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Stum & 3774 & 3000 & 9.24E-05 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Aldh1a1 & 5190 & 4006 & 1.42E-04 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Etv1 & 2097 & 1582 & 1.55E-04 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Syt13 & 5708 & 4378 & 1.57E-04 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Zfp365 & 5896 & 4907 & 1.69E-04 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Thrsp & 1440 & 1210 & 1.74E-04 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & 2310058D17Rik & 534 & 505 & 5.21E-04 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Myo5b & 1480 & 1160 & 5.44E-04 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Nudt4 & 8909 & 6478 & 6.03E-04 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Sphkap & 2039 & 1653 & 8.62E-04 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Synpr & 1001 & 831 & 1.14E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Cdo1 & 976 & 810 & 1.32E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Rps20 & 36917 & 23266 & 1.34E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Resp18 & 4095 & 3182 & 1.69E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Trim9 & 1997 & 1646 & 1.71E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Prkcq & 1487 & 1014 & 2.65E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Tunar & 1083 & 916 & 2.93E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Pik3r4 & 1370 & 1106 & 3.27E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Msh2 & 1565 & 1236 & 3.87E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Efhd2 & 4295 & 3551 & 4.35E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Slc35c2 & 902 & 785 & 6.24E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Golm1 & 665 & 565 & 6.83E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Igfbp4 & 2781 & 1787 & 7.46E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Rtn4 & 11088 & 8257 & 7.46E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Sema5a & 1319 & 1059 & 7.46E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Rgs4 & 6248 & 4578 & 7.54E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Ripor2 & 842 & 724 & 8.58E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Itpr1 & 3409 & 2816 & 8.70E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Ppp3ca & 38830 & 13587 & 8.70E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Elp6 & 2587 & 2431 & 9.49E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Klc1 & 1864 & 1553 & 9.72E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Gng5 & 1020 & 958 & 9.73E-03 & FALSE \\ 
      Hippocampus & Mouse 3 & Slide-seq & Lmo7 & 2167 & 1661 & 9.73E-03 & FALSE \\ 
      Cerebellum & Mouse 4 & Visium & Rps8 & 5145 & 2173 & 2.14E-03 & FALSE \\ 
      Cerebellum & Mouse 4 & Visium & Meg3 & 5375 & 2175 & 5.57E-03 & FALSE \\ 
       \hline
    \end{tabular}
    \end{table}

# Supp table ct spatial

``` r
st4 <-  hippo3_sp_ct |>
  mutate(tissue = 'Hippocampus', sample = 'Mouse 3', technology = 'Slide-seq') |>
  select(tissue, sample, technology, gene, cell_type,p_val,  xchr)  |>
  bind_rows(
    cere3_sp_ct |>
      mutate(tissue = 'Cerebellum', sample = 'Mouse 3', technology = 'Slide-seq') |>
      select(tissue, sample, technology, gene, cell_type,p_val,  xchr) 
  ) 
write.csv(st4, file = 'tables/04_supp_table_spatial_celltype_bias.csv')
```
