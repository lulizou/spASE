Plotting X-chromosome genes to see if any escape XCI
================

- <a href="#hippo-1" id="toc-hippo-1">Hippo 1</a>
- <a href="#hippo-2" id="toc-hippo-2">Hippo 2</a>
- <a href="#hippo-3" id="toc-hippo-3">Hippo 3</a>
- <a href="#cerebellum-3" id="toc-cerebellum-3">Cerebellum 3</a>
- <a href="#cerebellum-4" id="toc-cerebellum-4">Cerebellum 4</a>

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
library(Matrix)
library(data.table)
```


    Attaching package: 'data.table'

    The following objects are masked from 'package:dplyr':

        between, first, last

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

    The following objects are masked from 'package:data.table':

        first, second

    The following objects are masked from 'package:Matrix':

        expand, unname

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
library(tibble)
library(ggplot2)
library(ggthemes)
```

``` r
# Read in gencode to grab xchr genes
gencode <- import('results/gencode.vM10.annotation.gff3.gz')
xchr_genes <- unique(gencode$gene_name[which(seqnames(gencode)=='chrX')])
xchr_genes <- c(xchr_genes, 'Bex3')
```

## Hippo 1

``` r
hippo1 <- readRDS('results/rctd_hippo_1.rds')
maternal_counts_matrix_hippo <- hippo1@originalSpatialRNA@maternalCounts
paternal_counts_matrix_hippo <- hippo1@originalSpatialRNA@paternalCounts
coords_hippo <- hippo1@originalSpatialRNA@coords |> rownames_to_column()
```

``` r
xgenes <- rownames(maternal_counts_matrix_hippo)[(rownames(maternal_counts_matrix_hippo) %in% xchr_genes) & ((rowSums(maternal_counts_matrix_hippo) + rowSums(paternal_counts_matrix_hippo)) > 1000) & (rowSums(maternal_counts_matrix_hippo)>0 & rowSums(paternal_counts_matrix_hippo)>0)]
```

``` r
myfit_hippo <- spase(maternal_counts_matrix_hippo, paternal_counts_matrix_hippo, coords_hippo,cores=1,verbose=T, df=15, genes = xgenes)
```

    2057 genes pass min threshold of 100 pixels, 500 UMI

    running on 25 user-specified genes

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates


      |                                                                            
      |                                                                      |   0%
      |                                                                            
      |===                                                                   |   4%
      |                                                                            
      |======                                                                |   8%
      |                                                                            
      |========                                                              |  12%
      |                                                                            
      |===========                                                           |  16%
      |                                                                            
      |==============                                                        |  20%
      |                                                                            
      |=================                                                     |  24%
      |                                                                            
      |====================                                                  |  28%
      |                                                                            
      |======================                                                |  32%
      |                                                                            
      |=========================                                             |  36%
      |                                                                            
      |============================                                          |  40%
      |                                                                            
      |===============================                                       |  44%
      |                                                                            
      |==================================                                    |  48%
      |                                                                            
      |====================================                                  |  52%
      |                                                                            
      |=======================================                               |  56%
      |                                                                            
      |==========================================                            |  60%
      |                                                                            
      |=============================================                         |  64%
      |                                                                            
      |================================================                      |  68%
      |                                                                            
      |==================================================                    |  72%
      |                                                                            
      |=====================================================                 |  76%
      |                                                                            
      |========================================================              |  80%
      |                                                                            
      |===========================================================           |  84%
      |                                                                            
      |==============================================================        |  88%
      |                                                                            
      |================================================================      |  92%
      |                                                                            
      |===================================================================   |  96%
      |                                                                            
      |======================================================================| 100%

``` r
plotSpase(
  matrix1 = maternal_counts_matrix_hippo, 
  matrix2 = paternal_counts_matrix_hippo, 
  covariates = coords_hippo, 
  spasefit = myfit_hippo, 
  coords = coords_hippo |> select(x,y) |> sample_n(10000), 
  genes = xgenes,
  crosshairs = T,
  crosshairs_diag = F,
  point.size = 0.75,
  size.scale = F,
  theme = 'void',
  void = T,
  save = 'figures/08_xchr_hippo1'
)
```

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates

    [1] "Arhgef9"

    saved figures/08_xchr_hippo1-raw-Arhgef9.png

    saved figures/08_xchr_hippo1-smooth-Arhgef9.png

    saved figures/08_xchr_hippo1-zscore-Arhgef9.png

    saved figures/08_xchr_hippo1-x1-Arhgef9.png

    saved figures/08_xchr_hippo1-x2-Arhgef9.png

    [1] "Atrx"

    saved figures/08_xchr_hippo1-raw-Atrx.png

    saved figures/08_xchr_hippo1-smooth-Atrx.png

    saved figures/08_xchr_hippo1-zscore-Atrx.png

    saved figures/08_xchr_hippo1-x1-Atrx.png

    saved figures/08_xchr_hippo1-x2-Atrx.png

    [1] "Bex2"

    saved figures/08_xchr_hippo1-raw-Bex2.png

    saved figures/08_xchr_hippo1-smooth-Bex2.png

    saved figures/08_xchr_hippo1-zscore-Bex2.png

    saved figures/08_xchr_hippo1-x1-Bex2.png

    saved figures/08_xchr_hippo1-x2-Bex2.png

    [1] "Cdr1"

    saved figures/08_xchr_hippo1-raw-Cdr1.png

    saved figures/08_xchr_hippo1-smooth-Cdr1.png

    saved figures/08_xchr_hippo1-zscore-Cdr1.png

    saved figures/08_xchr_hippo1-x1-Cdr1.png

    saved figures/08_xchr_hippo1-x2-Cdr1.png

    [1] "Gnl3l"

    saved figures/08_xchr_hippo1-raw-Gnl3l.png

    saved figures/08_xchr_hippo1-smooth-Gnl3l.png

    saved figures/08_xchr_hippo1-zscore-Gnl3l.png

    saved figures/08_xchr_hippo1-x1-Gnl3l.png

    saved figures/08_xchr_hippo1-x2-Gnl3l.png

    [1] "Gpm6b"

    saved figures/08_xchr_hippo1-raw-Gpm6b.png

    saved figures/08_xchr_hippo1-smooth-Gpm6b.png

    saved figures/08_xchr_hippo1-zscore-Gpm6b.png

    saved figures/08_xchr_hippo1-x1-Gpm6b.png

    saved figures/08_xchr_hippo1-x2-Gpm6b.png

    [1] "Gprasp1"

    saved figures/08_xchr_hippo1-raw-Gprasp1.png

    saved figures/08_xchr_hippo1-smooth-Gprasp1.png

    saved figures/08_xchr_hippo1-zscore-Gprasp1.png

    saved figures/08_xchr_hippo1-x1-Gprasp1.png

    saved figures/08_xchr_hippo1-x2-Gprasp1.png

    [1] "Hprt"

    saved figures/08_xchr_hippo1-raw-Hprt.png

    saved figures/08_xchr_hippo1-smooth-Hprt.png

    saved figures/08_xchr_hippo1-zscore-Hprt.png

    saved figures/08_xchr_hippo1-x1-Hprt.png

    saved figures/08_xchr_hippo1-x2-Hprt.png

    [1] "Idh3g"

    saved figures/08_xchr_hippo1-raw-Idh3g.png

    saved figures/08_xchr_hippo1-smooth-Idh3g.png

    saved figures/08_xchr_hippo1-zscore-Idh3g.png

    saved figures/08_xchr_hippo1-x1-Idh3g.png

    saved figures/08_xchr_hippo1-x2-Idh3g.png

    [1] "Ids"

    saved figures/08_xchr_hippo1-raw-Ids.png

    saved figures/08_xchr_hippo1-smooth-Ids.png

    saved figures/08_xchr_hippo1-zscore-Ids.png

    saved figures/08_xchr_hippo1-x1-Ids.png

    saved figures/08_xchr_hippo1-x2-Ids.png

    [1] "Maged1"

    saved figures/08_xchr_hippo1-raw-Maged1.png

    saved figures/08_xchr_hippo1-smooth-Maged1.png

    saved figures/08_xchr_hippo1-zscore-Maged1.png

    saved figures/08_xchr_hippo1-x1-Maged1.png

    saved figures/08_xchr_hippo1-x2-Maged1.png

    [1] "Magee1"

    saved figures/08_xchr_hippo1-raw-Magee1.png

    saved figures/08_xchr_hippo1-smooth-Magee1.png

    saved figures/08_xchr_hippo1-zscore-Magee1.png

    saved figures/08_xchr_hippo1-x1-Magee1.png

    saved figures/08_xchr_hippo1-x2-Magee1.png

    [1] "Map7d2"

    saved figures/08_xchr_hippo1-raw-Map7d2.png

    saved figures/08_xchr_hippo1-smooth-Map7d2.png

    saved figures/08_xchr_hippo1-zscore-Map7d2.png

    saved figures/08_xchr_hippo1-x1-Map7d2.png

    saved figures/08_xchr_hippo1-x2-Map7d2.png

    [1] "Morf4l2"

    saved figures/08_xchr_hippo1-raw-Morf4l2.png

    saved figures/08_xchr_hippo1-smooth-Morf4l2.png

    saved figures/08_xchr_hippo1-zscore-Morf4l2.png

    saved figures/08_xchr_hippo1-x1-Morf4l2.png

    saved figures/08_xchr_hippo1-x2-Morf4l2.png

    [1] "Ndufb11"

    saved figures/08_xchr_hippo1-raw-Ndufb11.png

    saved figures/08_xchr_hippo1-smooth-Ndufb11.png

    saved figures/08_xchr_hippo1-zscore-Ndufb11.png

    saved figures/08_xchr_hippo1-x1-Ndufb11.png

    saved figures/08_xchr_hippo1-x2-Ndufb11.png

    [1] "Pcsk1n"

    saved figures/08_xchr_hippo1-raw-Pcsk1n.png

    saved figures/08_xchr_hippo1-smooth-Pcsk1n.png

    saved figures/08_xchr_hippo1-zscore-Pcsk1n.png

    saved figures/08_xchr_hippo1-x1-Pcsk1n.png

    saved figures/08_xchr_hippo1-x2-Pcsk1n.png

    [1] "Plp1"

    saved figures/08_xchr_hippo1-raw-Plp1.png

    saved figures/08_xchr_hippo1-smooth-Plp1.png

    saved figures/08_xchr_hippo1-zscore-Plp1.png

    saved figures/08_xchr_hippo1-x1-Plp1.png

    saved figures/08_xchr_hippo1-x2-Plp1.png

    [1] "Rps4x"

    saved figures/08_xchr_hippo1-raw-Rps4x.png

    saved figures/08_xchr_hippo1-smooth-Rps4x.png

    saved figures/08_xchr_hippo1-zscore-Rps4x.png

    saved figures/08_xchr_hippo1-x1-Rps4x.png

    saved figures/08_xchr_hippo1-x2-Rps4x.png

    [1] "Tceal3"

    saved figures/08_xchr_hippo1-raw-Tceal3.png

    saved figures/08_xchr_hippo1-smooth-Tceal3.png

    saved figures/08_xchr_hippo1-zscore-Tceal3.png

    saved figures/08_xchr_hippo1-x1-Tceal3.png

    saved figures/08_xchr_hippo1-x2-Tceal3.png

    [1] "Tceal5"

    saved figures/08_xchr_hippo1-raw-Tceal5.png

    saved figures/08_xchr_hippo1-smooth-Tceal5.png

    saved figures/08_xchr_hippo1-zscore-Tceal5.png

    saved figures/08_xchr_hippo1-x1-Tceal5.png

    saved figures/08_xchr_hippo1-x2-Tceal5.png

    [1] "Tceal6"

    saved figures/08_xchr_hippo1-raw-Tceal6.png

    saved figures/08_xchr_hippo1-smooth-Tceal6.png

    saved figures/08_xchr_hippo1-zscore-Tceal6.png

    saved figures/08_xchr_hippo1-x1-Tceal6.png

    saved figures/08_xchr_hippo1-x2-Tceal6.png

    [1] "Tspan7"

    saved figures/08_xchr_hippo1-raw-Tspan7.png

    saved figures/08_xchr_hippo1-smooth-Tspan7.png

    saved figures/08_xchr_hippo1-zscore-Tspan7.png

    saved figures/08_xchr_hippo1-x1-Tspan7.png

    saved figures/08_xchr_hippo1-x2-Tspan7.png

    [1] "Uba1"

    saved figures/08_xchr_hippo1-raw-Uba1.png

    saved figures/08_xchr_hippo1-smooth-Uba1.png

    saved figures/08_xchr_hippo1-zscore-Uba1.png

    saved figures/08_xchr_hippo1-x1-Uba1.png

    saved figures/08_xchr_hippo1-x2-Uba1.png

    [1] "Usp11"

    saved figures/08_xchr_hippo1-raw-Usp11.png

    saved figures/08_xchr_hippo1-smooth-Usp11.png

    saved figures/08_xchr_hippo1-zscore-Usp11.png

    saved figures/08_xchr_hippo1-x1-Usp11.png

    saved figures/08_xchr_hippo1-x2-Usp11.png

    [1] "Wdr13"

    saved figures/08_xchr_hippo1-raw-Wdr13.png

    saved figures/08_xchr_hippo1-smooth-Wdr13.png

    saved figures/08_xchr_hippo1-zscore-Wdr13.png

    saved figures/08_xchr_hippo1-x1-Wdr13.png

    saved figures/08_xchr_hippo1-x2-Wdr13.png

## Hippo 2

``` r
hippo2 <- readRDS('results/rctd_hippo_2.rds')
maternal_counts_matrix_hippo <- hippo2@originalSpatialRNA@maternalCounts
paternal_counts_matrix_hippo <- hippo2@originalSpatialRNA@paternalCounts
coords_hippo <- hippo2@originalSpatialRNA@coords |> rownames_to_column()
```

``` r
xgenes <- rownames(maternal_counts_matrix_hippo)[(rownames(maternal_counts_matrix_hippo) %in% xchr_genes) & ((rowSums(maternal_counts_matrix_hippo) + rowSums(paternal_counts_matrix_hippo)) > 1000) & (rowSums(maternal_counts_matrix_hippo)>0 & rowSums(paternal_counts_matrix_hippo)>0)]
xgenes <- c(xgenes, 'Morf4l2', 'Tceal5', 'Tceal6')
```

``` r
myfit_hippo <- spase(maternal_counts_matrix_hippo, paternal_counts_matrix_hippo, coords_hippo,cores=1,verbose=T, df=15, genes = xgenes, min.umi=100)
```

    3129 genes pass min threshold of 100 pixels, 100 UMI

    running on 8 user-specified genes

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates


      |                                                                            
      |                                                                      |   0%
      |                                                                            
      |=========                                                             |  12%
      |                                                                            
      |==================                                                    |  25%
      |                                                                            
      |==========================                                            |  38%
      |                                                                            
      |===================================                                   |  50%
      |                                                                            
      |============================================                          |  62%
      |                                                                            
      |====================================================                  |  75%
      |                                                                            
      |=============================================================         |  88%
      |                                                                            
      |======================================================================| 100%

``` r
plotSpase(
  matrix1 = maternal_counts_matrix_hippo, 
  matrix2 = paternal_counts_matrix_hippo, 
  covariates = coords_hippo, 
  spasefit = myfit_hippo, 
  coords = coords_hippo |> select(x,y) |> sample_n(10000), 
  genes = xgenes,
  crosshairs = T,
  crosshairs_diag = F,
  point.size = 0.75,
  size.scale = F,
  theme = 'void',
  void = T,
  save = 'figures/08_xchr_hippo2'
)
```

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates

    [1] "Ndufb11"

    saved figures/08_xchr_hippo2-raw-Ndufb11.png

    saved figures/08_xchr_hippo2-smooth-Ndufb11.png

    saved figures/08_xchr_hippo2-zscore-Ndufb11.png

    saved figures/08_xchr_hippo2-x1-Ndufb11.png

    saved figures/08_xchr_hippo2-x2-Ndufb11.png

    [1] "Plp1"

    saved figures/08_xchr_hippo2-raw-Plp1.png

    saved figures/08_xchr_hippo2-smooth-Plp1.png

    saved figures/08_xchr_hippo2-zscore-Plp1.png

    saved figures/08_xchr_hippo2-x1-Plp1.png

    saved figures/08_xchr_hippo2-x2-Plp1.png

    [1] "Syp"

    saved figures/08_xchr_hippo2-raw-Syp.png

    saved figures/08_xchr_hippo2-smooth-Syp.png

    saved figures/08_xchr_hippo2-zscore-Syp.png

    saved figures/08_xchr_hippo2-x1-Syp.png

    saved figures/08_xchr_hippo2-x2-Syp.png

    [1] "Tspan7"

    saved figures/08_xchr_hippo2-raw-Tspan7.png

    saved figures/08_xchr_hippo2-smooth-Tspan7.png

    saved figures/08_xchr_hippo2-zscore-Tspan7.png

    saved figures/08_xchr_hippo2-x1-Tspan7.png

    saved figures/08_xchr_hippo2-x2-Tspan7.png

    [1] "Uba1"

    saved figures/08_xchr_hippo2-raw-Uba1.png

    saved figures/08_xchr_hippo2-smooth-Uba1.png

    saved figures/08_xchr_hippo2-zscore-Uba1.png

    saved figures/08_xchr_hippo2-x1-Uba1.png

    saved figures/08_xchr_hippo2-x2-Uba1.png

    [1] "Morf4l2"

    Warning in plotSpase(matrix1 = maternal_counts_matrix_hippo, matrix2 =
    paternal_counts_matrix_hippo, : Morf4l2 did not converge; try lowering degrees
    of freedom. skipping for now

    [1] "Tceal5"

    saved figures/08_xchr_hippo2-raw-Tceal5.png

    saved figures/08_xchr_hippo2-smooth-Tceal5.png

    saved figures/08_xchr_hippo2-zscore-Tceal5.png

    saved figures/08_xchr_hippo2-x1-Tceal5.png

    saved figures/08_xchr_hippo2-x2-Tceal5.png

    [1] "Tceal6"

    saved figures/08_xchr_hippo2-raw-Tceal6.png

    saved figures/08_xchr_hippo2-smooth-Tceal6.png

    saved figures/08_xchr_hippo2-zscore-Tceal6.png

    saved figures/08_xchr_hippo2-x1-Tceal6.png

    saved figures/08_xchr_hippo2-x2-Tceal6.png

Smooth Morf4l2 with lower df to be able to compare to hippo1

``` r
myfit_hippo <- spase(maternal_counts_matrix_hippo, paternal_counts_matrix_hippo, coords_hippo,cores=1,verbose=T, df=5, genes = 'Morf4l2', min.umi=100)
```

    3129 genes pass min threshold of 100 pixels, 100 UMI

    running on 1 user-specified genes

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates


      |                                                                            
      |                                                                      |   0%
      |                                                                            
      |======================================================================| 100%

``` r
plotSpase(
  matrix1 = maternal_counts_matrix_hippo, 
  matrix2 = paternal_counts_matrix_hippo, 
  covariates = coords_hippo, 
  spasefit = myfit_hippo, 
  coords = coords_hippo |> select(x,y) |> sample_n(10000), 
  genes = 'Morf4l2',
  crosshairs = T,
  crosshairs_diag = F,
  point.size = 0.75,
  size.scale = F,
  theme = 'void',
  void = T,
  save = 'figures/08_xchr_hippo2'
)
```

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates

    [1] "Morf4l2"

    saved figures/08_xchr_hippo2-raw-Morf4l2.png

    saved figures/08_xchr_hippo2-smooth-Morf4l2.png

    saved figures/08_xchr_hippo2-zscore-Morf4l2.png

    saved figures/08_xchr_hippo2-x1-Morf4l2.png

    saved figures/08_xchr_hippo2-x2-Morf4l2.png

## Hippo 3

``` r
hippo3 <- readRDS('results/rctd_hippo_3.rds')
maternal_counts_matrix_hippo <- hippo3@originalSpatialRNA@maternalCounts
paternal_counts_matrix_hippo <- hippo3@originalSpatialRNA@paternalCounts
coords_hippo <- hippo3@originalSpatialRNA@coords |> rownames_to_column()
```

``` r
xgenes <- rownames(maternal_counts_matrix_hippo)[(rownames(maternal_counts_matrix_hippo) %in% xchr_genes) & ((rowSums(maternal_counts_matrix_hippo) + rowSums(paternal_counts_matrix_hippo)) > 1000) & (rowSums(maternal_counts_matrix_hippo)>0 & rowSums(paternal_counts_matrix_hippo)>0)]
```

``` r
myfit_hippo <- spase(maternal_counts_matrix_hippo, paternal_counts_matrix_hippo, coords_hippo,cores=1,verbose=T, df=15, genes = c('Arhgef9', 'Bex3', 'Cdkl5', 'Sat1', 'Syp', 'Tceal5'), min.umi=100) # subset to the ones i'm actually going to include in the supplemental figure; the rest are pretty much all maternal
```

    8922 genes pass min threshold of 100 pixels, 100 UMI

    running on 6 user-specified genes

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates


      |                                                                            
      |                                                                      |   0%
      |                                                                            
      |============                                                          |  17%
      |                                                                            
      |=======================                                               |  33%
      |                                                                            
      |===================================                                   |  50%
      |                                                                            
      |===============================================                       |  67%
      |                                                                            
      |==========================================================            |  83%
      |                                                                            
      |======================================================================| 100%

``` r
set.seed(5)
common_coords <- coords_hippo |> select(x,y) |> sample_n(10000)
plotSpase(
  matrix1 = maternal_counts_matrix_hippo, 
  matrix2 = paternal_counts_matrix_hippo, 
  covariates = coords_hippo, 
  spasefit = myfit_hippo, 
  coords = common_coords, 
  genes = c('Arhgef9', 'Bex3', 'Cdkl5', 'Sat1', 'Syp', 'Tceal5'),
  crosshairs = T,
  crosshairs_diag = F,
  point.size = 0.75,
  size.scale = F,
  theme = 'void',
  void = T,
  save = 'figures/08_xchr_hippo3'
)
```

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates

    [1] "Arhgef9"

    saved figures/08_xchr_hippo3-raw-Arhgef9.png

    saved figures/08_xchr_hippo3-smooth-Arhgef9.png

    saved figures/08_xchr_hippo3-zscore-Arhgef9.png

    saved figures/08_xchr_hippo3-x1-Arhgef9.png

    saved figures/08_xchr_hippo3-x2-Arhgef9.png

    [1] "Bex3"

    saved figures/08_xchr_hippo3-raw-Bex3.png

    saved figures/08_xchr_hippo3-smooth-Bex3.png

    saved figures/08_xchr_hippo3-zscore-Bex3.png

    saved figures/08_xchr_hippo3-x1-Bex3.png

    saved figures/08_xchr_hippo3-x2-Bex3.png

    [1] "Cdkl5"

    saved figures/08_xchr_hippo3-raw-Cdkl5.png

    saved figures/08_xchr_hippo3-smooth-Cdkl5.png

    saved figures/08_xchr_hippo3-zscore-Cdkl5.png

    saved figures/08_xchr_hippo3-x1-Cdkl5.png

    saved figures/08_xchr_hippo3-x2-Cdkl5.png

    [1] "Sat1"

    saved figures/08_xchr_hippo3-raw-Sat1.png

    saved figures/08_xchr_hippo3-smooth-Sat1.png

    saved figures/08_xchr_hippo3-zscore-Sat1.png

    saved figures/08_xchr_hippo3-x1-Sat1.png

    saved figures/08_xchr_hippo3-x2-Sat1.png

    [1] "Syp"

    saved figures/08_xchr_hippo3-raw-Syp.png

    saved figures/08_xchr_hippo3-smooth-Syp.png

    saved figures/08_xchr_hippo3-zscore-Syp.png

    saved figures/08_xchr_hippo3-x1-Syp.png

    saved figures/08_xchr_hippo3-x2-Syp.png

    [1] "Tceal5"

    saved figures/08_xchr_hippo3-raw-Tceal5.png

    saved figures/08_xchr_hippo3-smooth-Tceal5.png

    saved figures/08_xchr_hippo3-zscore-Tceal5.png

    saved figures/08_xchr_hippo3-x1-Tceal5.png

    saved figures/08_xchr_hippo3-x2-Tceal5.png

Extra genes for comparison to other samples

``` r
myfit_hippo <- spase(maternal_counts_matrix_hippo, paternal_counts_matrix_hippo, coords_hippo,cores=1,verbose=T, df=5, genes = c('Tceal3','Tceal6','Morf4l2'), min.umi=100)
```

    8922 genes pass min threshold of 100 pixels, 100 UMI

    running on 3 user-specified genes

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates


      |                                                                            
      |                                                                      |   0%
      |                                                                            
      |=======================                                               |  33%
      |                                                                            
      |===============================================                       |  67%
      |                                                                            
      |======================================================================| 100%

``` r
plotSpase(
  matrix1 = maternal_counts_matrix_hippo, 
  matrix2 = paternal_counts_matrix_hippo, 
  covariates = coords_hippo, 
  spasefit = myfit_hippo, 
  coords = common_coords, 
  genes = c('Tceal3','Tceal6','Morf4l2'),
  crosshairs = T,
  crosshairs_diag = F,
  point.size = 0.75,
  size.scale = F,
  theme = 'void',
  void = T,
  save = 'figures/08_xchr_hippo3'
)
```

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates

    [1] "Tceal3"

    saved figures/08_xchr_hippo3-raw-Tceal3.png

    saved figures/08_xchr_hippo3-smooth-Tceal3.png

    saved figures/08_xchr_hippo3-zscore-Tceal3.png

    saved figures/08_xchr_hippo3-x1-Tceal3.png

    saved figures/08_xchr_hippo3-x2-Tceal3.png

    [1] "Tceal6"

    saved figures/08_xchr_hippo3-raw-Tceal6.png

    saved figures/08_xchr_hippo3-smooth-Tceal6.png

    saved figures/08_xchr_hippo3-zscore-Tceal6.png

    saved figures/08_xchr_hippo3-x1-Tceal6.png

    saved figures/08_xchr_hippo3-x2-Tceal6.png

    [1] "Morf4l2"

    saved figures/08_xchr_hippo3-raw-Morf4l2.png

    saved figures/08_xchr_hippo3-smooth-Morf4l2.png

    saved figures/08_xchr_hippo3-zscore-Morf4l2.png

    saved figures/08_xchr_hippo3-x1-Morf4l2.png

    saved figures/08_xchr_hippo3-x2-Morf4l2.png

## Cerebellum 3

``` r
cere3 <- readRDS('results/rctd_cere_3.rds')
maternal_counts_matrix_cere <- cere3@originalSpatialRNA@maternalCounts
paternal_counts_matrix_cere <- cere3@originalSpatialRNA@paternalCounts
coords_cere <- cere3@originalSpatialRNA@coords |> rownames_to_column()
```

``` r
xgenes <- rownames(maternal_counts_matrix_cere)[(rownames(maternal_counts_matrix_cere) %in% xchr_genes) & ((rowSums(maternal_counts_matrix_cere) + rowSums(paternal_counts_matrix_cere)) > 1000) & (rowSums(maternal_counts_matrix_cere)>0 & rowSums(paternal_counts_matrix_cere)>0)]
```

``` r
myfit_cere <- spase(maternal_counts_matrix_cere, paternal_counts_matrix_cere, coords_cere,cores=1,verbose=T, df=15, genes = xgenes)
```

    4506 genes pass min threshold of 100 pixels, 500 UMI

    running on 55 user-specified genes

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates


      |                                                                            
      |                                                                      |   0%
      |                                                                            
      |=                                                                     |   2%
      |                                                                            
      |===                                                                   |   4%
      |                                                                            
      |====                                                                  |   5%
      |                                                                            
      |=====                                                                 |   7%
      |                                                                            
      |======                                                                |   9%
      |                                                                            
      |========                                                              |  11%
      |                                                                            
      |=========                                                             |  13%
      |                                                                            
      |==========                                                            |  15%
      |                                                                            
      |===========                                                           |  16%
      |                                                                            
      |=============                                                         |  18%
      |                                                                            
      |==============                                                        |  20%
      |                                                                            
      |===============                                                       |  22%
      |                                                                            
      |=================                                                     |  24%
      |                                                                            
      |==================                                                    |  25%
      |                                                                            
      |===================                                                   |  27%
      |                                                                            
      |====================                                                  |  29%
      |                                                                            
      |======================                                                |  31%
      |                                                                            
      |=======================                                               |  33%
      |                                                                            
      |========================                                              |  35%
      |                                                                            
      |=========================                                             |  36%
      |                                                                            
      |===========================                                           |  38%
      |                                                                            
      |============================                                          |  40%
      |                                                                            
      |=============================                                         |  42%
      |                                                                            
      |===============================                                       |  44%
      |                                                                            
      |================================                                      |  45%
      |                                                                            
      |=================================                                     |  47%
      |                                                                            
      |==================================                                    |  49%
      |                                                                            
      |====================================                                  |  51%
      |                                                                            
      |=====================================                                 |  53%
      |                                                                            
      |======================================                                |  55%
      |                                                                            
      |=======================================                               |  56%
      |                                                                            
      |=========================================                             |  58%
      |                                                                            
      |==========================================                            |  60%
      |                                                                            
      |===========================================                           |  62%
      |                                                                            
      |=============================================                         |  64%
      |                                                                            
      |==============================================                        |  65%
      |                                                                            
      |===============================================                       |  67%
      |                                                                            
      |================================================                      |  69%
      |                                                                            
      |==================================================                    |  71%
      |                                                                            
      |===================================================                   |  73%
      |                                                                            
      |====================================================                  |  75%
      |                                                                            
      |=====================================================                 |  76%
      |                                                                            
      |=======================================================               |  78%
      |                                                                            
      |========================================================              |  80%
      |                                                                            
      |=========================================================             |  82%
      |                                                                            
      |===========================================================           |  84%
      |                                                                            
      |============================================================          |  85%
      |                                                                            
      |=============================================================         |  87%
      |                                                                            
      |==============================================================        |  89%
      |                                                                            
      |================================================================      |  91%
      |                                                                            
      |=================================================================     |  93%
      |                                                                            
      |==================================================================    |  95%
      |                                                                            
      |===================================================================   |  96%
      |                                                                            
      |===================================================================== |  98%
      |                                                                            
      |======================================================================| 100%

``` r
plotSpase(
  matrix1 = maternal_counts_matrix_cere, 
  matrix2 = paternal_counts_matrix_cere, 
  covariates = coords_cere, 
  spasefit = myfit_cere, 
  coords = coords_cere |> select(x,y) |> sample_n(10000), 
  genes = xgenes,
  crosshairs = T,
  crosshairs_diag = F,
  point.size = 0.75,
  size.scale = F,
  theme = 'void',
  void = T,
  save = 'figures/08_xchr_cere3'
)
```

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates

    [1] "Araf"

    saved figures/08_xchr_cere3-raw-Araf.png

    saved figures/08_xchr_cere3-smooth-Araf.png

    saved figures/08_xchr_cere3-zscore-Araf.png

    saved figures/08_xchr_cere3-x1-Araf.png

    saved figures/08_xchr_cere3-x2-Araf.png

    [1] "Arhgef9"

    saved figures/08_xchr_cere3-raw-Arhgef9.png

    saved figures/08_xchr_cere3-smooth-Arhgef9.png

    saved figures/08_xchr_cere3-zscore-Arhgef9.png

    saved figures/08_xchr_cere3-x1-Arhgef9.png

    saved figures/08_xchr_cere3-x2-Arhgef9.png

    [1] "Armcx1"

    saved figures/08_xchr_cere3-raw-Armcx1.png

    saved figures/08_xchr_cere3-smooth-Armcx1.png

    saved figures/08_xchr_cere3-zscore-Armcx1.png

    saved figures/08_xchr_cere3-x1-Armcx1.png

    saved figures/08_xchr_cere3-x2-Armcx1.png

    [1] "Armcx3"

    saved figures/08_xchr_cere3-raw-Armcx3.png

    saved figures/08_xchr_cere3-smooth-Armcx3.png

    saved figures/08_xchr_cere3-zscore-Armcx3.png

    saved figures/08_xchr_cere3-x1-Armcx3.png

    saved figures/08_xchr_cere3-x2-Armcx3.png

    [1] "Arxes1"

    saved figures/08_xchr_cere3-raw-Arxes1.png

    saved figures/08_xchr_cere3-smooth-Arxes1.png

    saved figures/08_xchr_cere3-zscore-Arxes1.png

    saved figures/08_xchr_cere3-x1-Arxes1.png

    saved figures/08_xchr_cere3-x2-Arxes1.png

    [1] "Atrx"

    saved figures/08_xchr_cere3-raw-Atrx.png

    saved figures/08_xchr_cere3-smooth-Atrx.png

    saved figures/08_xchr_cere3-zscore-Atrx.png

    saved figures/08_xchr_cere3-x1-Atrx.png

    saved figures/08_xchr_cere3-x2-Atrx.png

    [1] "Bex2"

    saved figures/08_xchr_cere3-raw-Bex2.png

    saved figures/08_xchr_cere3-smooth-Bex2.png

    saved figures/08_xchr_cere3-zscore-Bex2.png

    saved figures/08_xchr_cere3-x1-Bex2.png

    saved figures/08_xchr_cere3-x2-Bex2.png

    [1] "Bex3"

    saved figures/08_xchr_cere3-raw-Bex3.png

    saved figures/08_xchr_cere3-smooth-Bex3.png

    saved figures/08_xchr_cere3-zscore-Bex3.png

    saved figures/08_xchr_cere3-x1-Bex3.png

    saved figures/08_xchr_cere3-x2-Bex3.png

    [1] "Cdk16"

    saved figures/08_xchr_cere3-raw-Cdk16.png

    saved figures/08_xchr_cere3-smooth-Cdk16.png

    saved figures/08_xchr_cere3-zscore-Cdk16.png

    saved figures/08_xchr_cere3-x1-Cdk16.png

    saved figures/08_xchr_cere3-x2-Cdk16.png

    [1] "Cdr1"

    saved figures/08_xchr_cere3-raw-Cdr1.png

    saved figures/08_xchr_cere3-smooth-Cdr1.png

    saved figures/08_xchr_cere3-zscore-Cdr1.png

    saved figures/08_xchr_cere3-x1-Cdr1.png

    saved figures/08_xchr_cere3-x2-Cdr1.png

    [1] "Cnksr2"

    saved figures/08_xchr_cere3-raw-Cnksr2.png

    saved figures/08_xchr_cere3-smooth-Cnksr2.png

    saved figures/08_xchr_cere3-zscore-Cnksr2.png

    saved figures/08_xchr_cere3-x1-Cnksr2.png

    saved figures/08_xchr_cere3-x2-Cnksr2.png

    [1] "Ctps2"

    saved figures/08_xchr_cere3-raw-Ctps2.png

    saved figures/08_xchr_cere3-smooth-Ctps2.png

    saved figures/08_xchr_cere3-zscore-Ctps2.png

    saved figures/08_xchr_cere3-x1-Ctps2.png

    saved figures/08_xchr_cere3-x2-Ctps2.png

    [1] "Dkc1"

    saved figures/08_xchr_cere3-raw-Dkc1.png

    saved figures/08_xchr_cere3-smooth-Dkc1.png

    saved figures/08_xchr_cere3-zscore-Dkc1.png

    saved figures/08_xchr_cere3-x1-Dkc1.png

    saved figures/08_xchr_cere3-x2-Dkc1.png

    [1] "Eif1ax"

    saved figures/08_xchr_cere3-raw-Eif1ax.png

    saved figures/08_xchr_cere3-smooth-Eif1ax.png

    saved figures/08_xchr_cere3-zscore-Eif1ax.png

    saved figures/08_xchr_cere3-x1-Eif1ax.png

    saved figures/08_xchr_cere3-x2-Eif1ax.png

    [1] "Fundc2"

    saved figures/08_xchr_cere3-raw-Fundc2.png

    saved figures/08_xchr_cere3-smooth-Fundc2.png

    saved figures/08_xchr_cere3-zscore-Fundc2.png

    saved figures/08_xchr_cere3-x1-Fundc2.png

    saved figures/08_xchr_cere3-x2-Fundc2.png

    [1] "Gnl3l"

    saved figures/08_xchr_cere3-raw-Gnl3l.png

    saved figures/08_xchr_cere3-smooth-Gnl3l.png

    saved figures/08_xchr_cere3-zscore-Gnl3l.png

    saved figures/08_xchr_cere3-x1-Gnl3l.png

    saved figures/08_xchr_cere3-x2-Gnl3l.png

    [1] "Gpm6b"

    saved figures/08_xchr_cere3-raw-Gpm6b.png

    saved figures/08_xchr_cere3-smooth-Gpm6b.png

    saved figures/08_xchr_cere3-zscore-Gpm6b.png

    saved figures/08_xchr_cere3-x1-Gpm6b.png

    saved figures/08_xchr_cere3-x2-Gpm6b.png

    [1] "Gprasp1"

    saved figures/08_xchr_cere3-raw-Gprasp1.png

    saved figures/08_xchr_cere3-smooth-Gprasp1.png

    saved figures/08_xchr_cere3-zscore-Gprasp1.png

    saved figures/08_xchr_cere3-x1-Gprasp1.png

    saved figures/08_xchr_cere3-x2-Gprasp1.png

    [1] "Hprt"

    saved figures/08_xchr_cere3-raw-Hprt.png

    saved figures/08_xchr_cere3-smooth-Hprt.png

    saved figures/08_xchr_cere3-zscore-Hprt.png

    saved figures/08_xchr_cere3-x1-Hprt.png

    saved figures/08_xchr_cere3-x2-Hprt.png

    [1] "Hsd17b10"

    saved figures/08_xchr_cere3-raw-Hsd17b10.png

    saved figures/08_xchr_cere3-smooth-Hsd17b10.png

    saved figures/08_xchr_cere3-zscore-Hsd17b10.png

    saved figures/08_xchr_cere3-x1-Hsd17b10.png

    saved figures/08_xchr_cere3-x2-Hsd17b10.png

    [1] "Htatsf1"

    saved figures/08_xchr_cere3-raw-Htatsf1.png

    saved figures/08_xchr_cere3-smooth-Htatsf1.png

    saved figures/08_xchr_cere3-zscore-Htatsf1.png

    saved figures/08_xchr_cere3-x1-Htatsf1.png

    saved figures/08_xchr_cere3-x2-Htatsf1.png

    [1] "Idh3g"

    saved figures/08_xchr_cere3-raw-Idh3g.png

    saved figures/08_xchr_cere3-smooth-Idh3g.png

    saved figures/08_xchr_cere3-zscore-Idh3g.png

    saved figures/08_xchr_cere3-x1-Idh3g.png

    saved figures/08_xchr_cere3-x2-Idh3g.png

    [1] "Ids"

    saved figures/08_xchr_cere3-raw-Ids.png

    saved figures/08_xchr_cere3-smooth-Ids.png

    saved figures/08_xchr_cere3-zscore-Ids.png

    saved figures/08_xchr_cere3-x1-Ids.png

    saved figures/08_xchr_cere3-x2-Ids.png

    [1] "Irak1"

    saved figures/08_xchr_cere3-raw-Irak1.png

    saved figures/08_xchr_cere3-smooth-Irak1.png

    saved figures/08_xchr_cere3-zscore-Irak1.png

    saved figures/08_xchr_cere3-x1-Irak1.png

    saved figures/08_xchr_cere3-x2-Irak1.png

    [1] "Maged1"

    saved figures/08_xchr_cere3-raw-Maged1.png

    saved figures/08_xchr_cere3-smooth-Maged1.png

    saved figures/08_xchr_cere3-zscore-Maged1.png

    saved figures/08_xchr_cere3-x1-Maged1.png

    saved figures/08_xchr_cere3-x2-Maged1.png

    [1] "Mageh1"

    saved figures/08_xchr_cere3-raw-Mageh1.png

    saved figures/08_xchr_cere3-smooth-Mageh1.png

    saved figures/08_xchr_cere3-zscore-Mageh1.png

    saved figures/08_xchr_cere3-x1-Mageh1.png

    saved figures/08_xchr_cere3-x2-Mageh1.png

    [1] "Map7d2"

    saved figures/08_xchr_cere3-raw-Map7d2.png

    saved figures/08_xchr_cere3-smooth-Map7d2.png

    saved figures/08_xchr_cere3-zscore-Map7d2.png

    saved figures/08_xchr_cere3-x1-Map7d2.png

    saved figures/08_xchr_cere3-x2-Map7d2.png

    [1] "Mcts1"

    saved figures/08_xchr_cere3-raw-Mcts1.png

    saved figures/08_xchr_cere3-smooth-Mcts1.png

    saved figures/08_xchr_cere3-zscore-Mcts1.png

    saved figures/08_xchr_cere3-x1-Mcts1.png

    saved figures/08_xchr_cere3-x2-Mcts1.png

    [1] "Ndufb11"

    saved figures/08_xchr_cere3-raw-Ndufb11.png

    saved figures/08_xchr_cere3-smooth-Ndufb11.png

    saved figures/08_xchr_cere3-zscore-Ndufb11.png

    saved figures/08_xchr_cere3-x1-Ndufb11.png

    saved figures/08_xchr_cere3-x2-Ndufb11.png

    [1] "Ogt"

    saved figures/08_xchr_cere3-raw-Ogt.png

    saved figures/08_xchr_cere3-smooth-Ogt.png

    saved figures/08_xchr_cere3-zscore-Ogt.png

    saved figures/08_xchr_cere3-x1-Ogt.png

    saved figures/08_xchr_cere3-x2-Ogt.png

    [1] "Pcsk1n"

    saved figures/08_xchr_cere3-raw-Pcsk1n.png

    saved figures/08_xchr_cere3-smooth-Pcsk1n.png

    saved figures/08_xchr_cere3-zscore-Pcsk1n.png

    saved figures/08_xchr_cere3-x1-Pcsk1n.png

    saved figures/08_xchr_cere3-x2-Pcsk1n.png

    [1] "Pdzd4"

    saved figures/08_xchr_cere3-raw-Pdzd4.png

    saved figures/08_xchr_cere3-smooth-Pdzd4.png

    saved figures/08_xchr_cere3-zscore-Pdzd4.png

    saved figures/08_xchr_cere3-x1-Pdzd4.png

    saved figures/08_xchr_cere3-x2-Pdzd4.png

    [1] "Pgrmc1"

    saved figures/08_xchr_cere3-raw-Pgrmc1.png

    saved figures/08_xchr_cere3-smooth-Pgrmc1.png

    saved figures/08_xchr_cere3-zscore-Pgrmc1.png

    saved figures/08_xchr_cere3-x1-Pgrmc1.png

    saved figures/08_xchr_cere3-x2-Pgrmc1.png

    [1] "Pja1"

    saved figures/08_xchr_cere3-raw-Pja1.png

    saved figures/08_xchr_cere3-smooth-Pja1.png

    saved figures/08_xchr_cere3-zscore-Pja1.png

    saved figures/08_xchr_cere3-x1-Pja1.png

    saved figures/08_xchr_cere3-x2-Pja1.png

    [1] "Plp1"

    saved figures/08_xchr_cere3-raw-Plp1.png

    saved figures/08_xchr_cere3-smooth-Plp1.png

    saved figures/08_xchr_cere3-zscore-Plp1.png

    saved figures/08_xchr_cere3-x1-Plp1.png

    saved figures/08_xchr_cere3-x2-Plp1.png

    [1] "Prps1"

    saved figures/08_xchr_cere3-raw-Prps1.png

    saved figures/08_xchr_cere3-smooth-Prps1.png

    saved figures/08_xchr_cere3-zscore-Prps1.png

    saved figures/08_xchr_cere3-x1-Prps1.png

    saved figures/08_xchr_cere3-x2-Prps1.png

    [1] "Rps4x"

    saved figures/08_xchr_cere3-raw-Rps4x.png

    saved figures/08_xchr_cere3-smooth-Rps4x.png

    saved figures/08_xchr_cere3-zscore-Rps4x.png

    saved figures/08_xchr_cere3-x1-Rps4x.png

    saved figures/08_xchr_cere3-x2-Rps4x.png

    [1] "Sat1"

    saved figures/08_xchr_cere3-raw-Sat1.png

    saved figures/08_xchr_cere3-smooth-Sat1.png

    saved figures/08_xchr_cere3-zscore-Sat1.png

    saved figures/08_xchr_cere3-x1-Sat1.png

    saved figures/08_xchr_cere3-x2-Sat1.png

    [1] "Sh3bgrl"

    saved figures/08_xchr_cere3-raw-Sh3bgrl.png

    saved figures/08_xchr_cere3-smooth-Sh3bgrl.png

    saved figures/08_xchr_cere3-zscore-Sh3bgrl.png

    saved figures/08_xchr_cere3-x1-Sh3bgrl.png

    saved figures/08_xchr_cere3-x2-Sh3bgrl.png

    [1] "Slc25a14"

    saved figures/08_xchr_cere3-raw-Slc25a14.png

    saved figures/08_xchr_cere3-smooth-Slc25a14.png

    saved figures/08_xchr_cere3-zscore-Slc25a14.png

    saved figures/08_xchr_cere3-x1-Slc25a14.png

    saved figures/08_xchr_cere3-x2-Slc25a14.png

    [1] "Slc6a8"

    saved figures/08_xchr_cere3-raw-Slc6a8.png

    saved figures/08_xchr_cere3-smooth-Slc6a8.png

    saved figures/08_xchr_cere3-zscore-Slc6a8.png

    saved figures/08_xchr_cere3-x1-Slc6a8.png

    saved figures/08_xchr_cere3-x2-Slc6a8.png

    [1] "Sms"

    saved figures/08_xchr_cere3-raw-Sms.png

    saved figures/08_xchr_cere3-smooth-Sms.png

    saved figures/08_xchr_cere3-zscore-Sms.png

    saved figures/08_xchr_cere3-x1-Sms.png

    saved figures/08_xchr_cere3-x2-Sms.png

    [1] "Syp"

    saved figures/08_xchr_cere3-raw-Syp.png

    saved figures/08_xchr_cere3-smooth-Syp.png

    saved figures/08_xchr_cere3-zscore-Syp.png

    saved figures/08_xchr_cere3-x1-Syp.png

    saved figures/08_xchr_cere3-x2-Syp.png

    [1] "Tceal3"

    saved figures/08_xchr_cere3-raw-Tceal3.png

    saved figures/08_xchr_cere3-smooth-Tceal3.png

    saved figures/08_xchr_cere3-zscore-Tceal3.png

    saved figures/08_xchr_cere3-x1-Tceal3.png

    saved figures/08_xchr_cere3-x2-Tceal3.png

    [1] "Tceal5"

    saved figures/08_xchr_cere3-raw-Tceal5.png

    saved figures/08_xchr_cere3-smooth-Tceal5.png

    saved figures/08_xchr_cere3-zscore-Tceal5.png

    saved figures/08_xchr_cere3-x1-Tceal5.png

    saved figures/08_xchr_cere3-x2-Tceal5.png

    [1] "Tceal6"

    saved figures/08_xchr_cere3-raw-Tceal6.png

    saved figures/08_xchr_cere3-smooth-Tceal6.png

    saved figures/08_xchr_cere3-zscore-Tceal6.png

    saved figures/08_xchr_cere3-x1-Tceal6.png

    saved figures/08_xchr_cere3-x2-Tceal6.png

    [1] "Tceal8"

    saved figures/08_xchr_cere3-raw-Tceal8.png

    saved figures/08_xchr_cere3-smooth-Tceal8.png

    saved figures/08_xchr_cere3-zscore-Tceal8.png

    saved figures/08_xchr_cere3-x1-Tceal8.png

    saved figures/08_xchr_cere3-x2-Tceal8.png

    [1] "Tmem47"

    saved figures/08_xchr_cere3-raw-Tmem47.png

    saved figures/08_xchr_cere3-smooth-Tmem47.png

    saved figures/08_xchr_cere3-zscore-Tmem47.png

    saved figures/08_xchr_cere3-x1-Tmem47.png

    saved figures/08_xchr_cere3-x2-Tmem47.png

    [1] "Trappc2"

    saved figures/08_xchr_cere3-raw-Trappc2.png

    saved figures/08_xchr_cere3-smooth-Trappc2.png

    saved figures/08_xchr_cere3-zscore-Trappc2.png

    saved figures/08_xchr_cere3-x1-Trappc2.png

    saved figures/08_xchr_cere3-x2-Trappc2.png

    [1] "Tspan7"

    saved figures/08_xchr_cere3-raw-Tspan7.png

    saved figures/08_xchr_cere3-smooth-Tspan7.png

    saved figures/08_xchr_cere3-zscore-Tspan7.png

    saved figures/08_xchr_cere3-x1-Tspan7.png

    saved figures/08_xchr_cere3-x2-Tspan7.png

    [1] "Tsr2"

    saved figures/08_xchr_cere3-raw-Tsr2.png

    saved figures/08_xchr_cere3-smooth-Tsr2.png

    saved figures/08_xchr_cere3-zscore-Tsr2.png

    saved figures/08_xchr_cere3-x1-Tsr2.png

    saved figures/08_xchr_cere3-x2-Tsr2.png

    [1] "Uba1"

    saved figures/08_xchr_cere3-raw-Uba1.png

    saved figures/08_xchr_cere3-smooth-Uba1.png

    saved figures/08_xchr_cere3-zscore-Uba1.png

    saved figures/08_xchr_cere3-x1-Uba1.png

    saved figures/08_xchr_cere3-x2-Uba1.png

    [1] "Vma21"

    saved figures/08_xchr_cere3-raw-Vma21.png

    saved figures/08_xchr_cere3-smooth-Vma21.png

    saved figures/08_xchr_cere3-zscore-Vma21.png

    saved figures/08_xchr_cere3-x1-Vma21.png

    saved figures/08_xchr_cere3-x2-Vma21.png

    [1] "Wdr13"

    saved figures/08_xchr_cere3-raw-Wdr13.png

    saved figures/08_xchr_cere3-smooth-Wdr13.png

    saved figures/08_xchr_cere3-zscore-Wdr13.png

    saved figures/08_xchr_cere3-x1-Wdr13.png

    saved figures/08_xchr_cere3-x2-Wdr13.png

    [1] "Xist"

    saved figures/08_xchr_cere3-raw-Xist.png

    saved figures/08_xchr_cere3-smooth-Xist.png

    saved figures/08_xchr_cere3-zscore-Xist.png

    saved figures/08_xchr_cere3-x1-Xist.png

    saved figures/08_xchr_cere3-x2-Xist.png

Gene I want to compare to other samples

``` r
myfit_cere <- spase(maternal_counts_matrix_cere, paternal_counts_matrix_cere, coords_cere,cores=1,verbose=T, df=5, genes = c('Morf4l2'), min.umi=100)
```

    8937 genes pass min threshold of 100 pixels, 100 UMI

    running on 1 user-specified genes

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates


      |                                                                            
      |                                                                      |   0%
      |                                                                            
      |======================================================================| 100%

``` r
plotSpase(
  matrix1 = maternal_counts_matrix_cere, 
  matrix2 = paternal_counts_matrix_cere, 
  covariates = coords_cere, 
  spasefit = myfit_cere, 
  coords = coords_cere |> select(x,y) |> sample_n(10000), 
  genes = 'Morf4l2',
  crosshairs = T,
  crosshairs_diag = F,
  point.size = 0.75,
  size.scale = F,
  theme = 'void',
  void = T,
  save = 'figures/08_xchr_cere3'
)
```

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates

    [1] "Morf4l2"

    saved figures/08_xchr_cere3-raw-Morf4l2.png

    saved figures/08_xchr_cere3-smooth-Morf4l2.png

    saved figures/08_xchr_cere3-zscore-Morf4l2.png

    saved figures/08_xchr_cere3-x1-Morf4l2.png

    saved figures/08_xchr_cere3-x2-Morf4l2.png

## Cerebellum 4

``` r
cere4 <- readRDS('results/rctd_cere_4_visium.rds')
maternal_counts_matrix_cere <- cere4@originalSpatialRNA@maternalCounts
paternal_counts_matrix_cere <- cere4@originalSpatialRNA@paternalCounts
coords_cere <- cere4@originalSpatialRNA@coords |> rownames_to_column()
```

``` r
xgenes <- rownames(maternal_counts_matrix_cere)[(rownames(maternal_counts_matrix_cere) %in% xchr_genes) & ((rowSums(maternal_counts_matrix_cere) + rowSums(paternal_counts_matrix_cere)) > 1000) & (rowSums(maternal_counts_matrix_cere)>0 & rowSums(paternal_counts_matrix_cere)>0)]
```

``` r
myfit_cere <- spase(maternal_counts_matrix_cere, paternal_counts_matrix_cere, coords_cere,cores=1,verbose=T, df=15, genes = xgenes)
```

    1321 genes pass min threshold of 100 pixels, 500 UMI

    running on 12 user-specified genes

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates


      |                                                                            
      |                                                                      |   0%
      |                                                                            
      |======                                                                |   8%
      |                                                                            
      |============                                                          |  17%
      |                                                                            
      |==================                                                    |  25%
      |                                                                            
      |=======================                                               |  33%
      |                                                                            
      |=============================                                         |  42%
      |                                                                            
      |===================================                                   |  50%
      |                                                                            
      |=========================================                             |  58%
      |                                                                            
      |===============================================                       |  67%
      |                                                                            
      |====================================================                  |  75%
      |                                                                            
      |==========================================================            |  83%
      |                                                                            
      |================================================================      |  92%
      |                                                                            
      |======================================================================| 100%

``` r
plotSpase(
  matrix1 = maternal_counts_matrix_cere, 
  matrix2 = paternal_counts_matrix_cere, 
  covariates = coords_cere, 
  spasefit = myfit_cere, 
  coords = coords_cere |> select(x,y) |> filter(y > 7500), # not a lot of reads at the bottom 
  genes = xgenes,
  crosshairs = T,
  crosshairs_diag = F,
  point.size = 0.75,
  size.scale = F,
  theme = 'void',
  void = T,
  save = 'figures/08_xchr_cere4'
)
```

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates

    [1] "Bex2"

    saved figures/08_xchr_cere4-raw-Bex2.png

    saved figures/08_xchr_cere4-smooth-Bex2.png

    saved figures/08_xchr_cere4-zscore-Bex2.png

    saved figures/08_xchr_cere4-x1-Bex2.png

    saved figures/08_xchr_cere4-x2-Bex2.png

    [1] "Cdr1"

    saved figures/08_xchr_cere4-raw-Cdr1.png

    saved figures/08_xchr_cere4-smooth-Cdr1.png

    saved figures/08_xchr_cere4-zscore-Cdr1.png

    saved figures/08_xchr_cere4-x1-Cdr1.png

    saved figures/08_xchr_cere4-x2-Cdr1.png

    [1] "Gprasp1"

    saved figures/08_xchr_cere4-raw-Gprasp1.png

    saved figures/08_xchr_cere4-smooth-Gprasp1.png

    saved figures/08_xchr_cere4-zscore-Gprasp1.png

    saved figures/08_xchr_cere4-x1-Gprasp1.png

    saved figures/08_xchr_cere4-x2-Gprasp1.png

    [1] "Idh3g"

    saved figures/08_xchr_cere4-raw-Idh3g.png

    saved figures/08_xchr_cere4-smooth-Idh3g.png

    saved figures/08_xchr_cere4-zscore-Idh3g.png

    saved figures/08_xchr_cere4-x1-Idh3g.png

    saved figures/08_xchr_cere4-x2-Idh3g.png

    [1] "Ids"

    saved figures/08_xchr_cere4-raw-Ids.png

    saved figures/08_xchr_cere4-smooth-Ids.png

    saved figures/08_xchr_cere4-zscore-Ids.png

    saved figures/08_xchr_cere4-x1-Ids.png

    saved figures/08_xchr_cere4-x2-Ids.png

    [1] "Maged1"

    saved figures/08_xchr_cere4-raw-Maged1.png

    saved figures/08_xchr_cere4-smooth-Maged1.png

    saved figures/08_xchr_cere4-zscore-Maged1.png

    saved figures/08_xchr_cere4-x1-Maged1.png

    saved figures/08_xchr_cere4-x2-Maged1.png

    [1] "Pcsk1n"

    saved figures/08_xchr_cere4-raw-Pcsk1n.png

    saved figures/08_xchr_cere4-smooth-Pcsk1n.png

    saved figures/08_xchr_cere4-zscore-Pcsk1n.png

    saved figures/08_xchr_cere4-x1-Pcsk1n.png

    saved figures/08_xchr_cere4-x2-Pcsk1n.png

    [1] "Plp1"

    saved figures/08_xchr_cere4-raw-Plp1.png

    saved figures/08_xchr_cere4-smooth-Plp1.png

    saved figures/08_xchr_cere4-zscore-Plp1.png

    saved figures/08_xchr_cere4-x1-Plp1.png

    saved figures/08_xchr_cere4-x2-Plp1.png

    [1] "Rps4x"

    saved figures/08_xchr_cere4-raw-Rps4x.png

    saved figures/08_xchr_cere4-smooth-Rps4x.png

    saved figures/08_xchr_cere4-zscore-Rps4x.png

    saved figures/08_xchr_cere4-x1-Rps4x.png

    saved figures/08_xchr_cere4-x2-Rps4x.png

    [1] "Tceal5"

    saved figures/08_xchr_cere4-raw-Tceal5.png

    saved figures/08_xchr_cere4-smooth-Tceal5.png

    saved figures/08_xchr_cere4-zscore-Tceal5.png

    saved figures/08_xchr_cere4-x1-Tceal5.png

    saved figures/08_xchr_cere4-x2-Tceal5.png

    [1] "Tspan7"

    saved figures/08_xchr_cere4-raw-Tspan7.png

    saved figures/08_xchr_cere4-smooth-Tspan7.png

    saved figures/08_xchr_cere4-zscore-Tspan7.png

    saved figures/08_xchr_cere4-x1-Tspan7.png

    saved figures/08_xchr_cere4-x2-Tspan7.png

    [1] "Uba1"

    saved figures/08_xchr_cere4-raw-Uba1.png

    saved figures/08_xchr_cere4-smooth-Uba1.png

    saved figures/08_xchr_cere4-zscore-Uba1.png

    saved figures/08_xchr_cere4-x1-Uba1.png

    saved figures/08_xchr_cere4-x2-Uba1.png

Add in ones to compare across samples; using lower df since fewer
spots/UMIs sampled for these genes

``` r
myfit_cere <- spase(maternal_counts_matrix_cere, paternal_counts_matrix_cere, coords_cere,cores=1,verbose=T, df=5, genes = c('Morf4l2', 'Tceal3', 'Tceal6', 'Xist'), min.umi=100)
```

    4848 genes pass min threshold of 100 pixels, 100 UMI

    running on 4 user-specified genes

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates


      |                                                                            
      |                                                                      |   0%
      |                                                                            
      |==================                                                    |  25%
      |                                                                            
      |===================================                                   |  50%
      |                                                                            
      |====================================================                  |  75%
      |                                                                            
      |======================================================================| 100%

``` r
plotSpase(
  matrix1 = maternal_counts_matrix_cere, 
  matrix2 = paternal_counts_matrix_cere, 
  covariates = coords_cere, 
  spasefit = myfit_cere, 
  coords = coords_cere |> select(x,y) |> filter(y > 7500), # not a lot of reads at the bottom 
  genes = c('Morf4l2', 'Tceal3', 'Tceal6','Xist'),
  crosshairs = T,
  crosshairs_diag = F,
  point.size = 0.75,
  size.scale = F,
  theme = 'void',
  void = T,
  save = 'figures/08_xchr_cere4'
)
```

    found 3 columns in covariates; going to assume that first column is pixel names, 2nd and 3rd column are 2D coordinates

    [1] "Morf4l2"

    saved figures/08_xchr_cere4-raw-Morf4l2.png

    saved figures/08_xchr_cere4-smooth-Morf4l2.png

    saved figures/08_xchr_cere4-zscore-Morf4l2.png

    saved figures/08_xchr_cere4-x1-Morf4l2.png

    saved figures/08_xchr_cere4-x2-Morf4l2.png

    [1] "Tceal3"

    saved figures/08_xchr_cere4-raw-Tceal3.png

    saved figures/08_xchr_cere4-smooth-Tceal3.png

    saved figures/08_xchr_cere4-zscore-Tceal3.png

    saved figures/08_xchr_cere4-x1-Tceal3.png

    saved figures/08_xchr_cere4-x2-Tceal3.png

    [1] "Tceal6"

    saved figures/08_xchr_cere4-raw-Tceal6.png

    saved figures/08_xchr_cere4-smooth-Tceal6.png

    saved figures/08_xchr_cere4-zscore-Tceal6.png

    saved figures/08_xchr_cere4-x1-Tceal6.png

    saved figures/08_xchr_cere4-x2-Tceal6.png

    [1] "Xist"

    saved figures/08_xchr_cere4-raw-Xist.png

    saved figures/08_xchr_cere4-smooth-Xist.png

    saved figures/08_xchr_cere4-zscore-Xist.png

    saved figures/08_xchr_cere4-x1-Xist.png

    saved figures/08_xchr_cere4-x2-Xist.png
