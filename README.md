spASE
================

- <a href="#installation" id="toc-installation">Installation</a>
- <a href="#quick-start" id="toc-quick-start">Quick start</a>
  - <a href="#hypothesis-testing" id="toc-hypothesis-testing">Hypothesis
    testing</a>
  - <a href="#visualization" id="toc-visualization">Visualization</a>
- <a href="#reproducibility" id="toc-reproducibility">Reproducibility</a>
- <a href="#citation" id="toc-citation">Citation</a>

R package for detecting 2D spatial patterns of allele-specific
expression in spatial transcriptomics data while controlling for cell
type. Provides functions for estimating, plotting, and testing for
significant ASE patterns. This repository also contains the code to
reproduce the results in the manuscript.

## Installation

``` r
devtools::install_github("lulizou/spASE")
```

## Quick start

The main input to spASE includes:

- maternal counts matrix (rows are genes, columns are spots or cells)
- paternal counts matrix
- spatial coordinates.

Note that spASE currently modifies many functions from the `spacexr`
package, so it is always good to load them in order so that the spASE
version is used:

``` r
library(spacexr)
library(spASE)
```

### Hypothesis testing

Testing for…

1.  Overall bias (i.e. pseudo-bulk):

``` r
bias_results <- scase(
  matrix1 = maternal_counts_matrix,
  matrix2 = paternal_counts_matrix,
  min.cells = 100, # set to something reasonable
  cores = 1, # can add more
  verbose = T # to enable progress bar
)
```

Note that this can also be run on a subset of genes by specifying them
as a string vector using the `genes = c(...)` argument.

2.  Overall spatial pattern (i.e. no cell type effect):

``` r
overall_spatial_results <- spase(
  matrix1 = maternal_counts_matrix,
  matrix2 = paternal_counts_matrix,
  covariates = coords,
  min.pixels = 100, # set to something reasonable
  cores = 1, # can add more
  verbose = T # to enable progress bar
)
```

Note that `coords` must have the first column as the pixel ID (matching
the column names of the count matrices) and the second two columns as
the spatial coordinates.

### Visualization

Once you have `overall_spatial_results`, we can visualize using the
`plotSpase` function as follows:

``` r
plotSpase(
  matrix1 = maternal_counts_matrix,
  matrix2 = paternal_counts_matrix,
  covariates = coords,
  spasefit = overall_spatial_results,
  coords = NULL, # to downsample could set e.g coords = coords |> sample_n(1e4)
)
```

## Reproducibility

The data is available in the Broad single cell portal at
[SCP1692](https://singlecell.broadinstitute.org/single_cell/study/SCP1692/detection-of-allele-specific-expression-in-spatial-transcriptomics-with-spase).

To reproduce the results in the manuscript, see the
[analysis](https://github.com/lulizou/spASE/tree/master/analysis)
folder. To see the knitted Quarto documents, see:

- [00_preprocess.md](https://github.com/lulizou/spASE/blob/master/analysis/00_preprocess.md) -
  details on the alignment and custom scripts for running RCTD, C-SIDE,
  and spASE
- [01_cell_type_maps.md](https://github.com/lulizou/spASE/blob/master/analysis/01_cell_type_maps.md) -
  cell type plots for each sample
- [02_compare_visium_slideseq_cerebellum.md](https://github.com/lulizou/spASE/blob/master/analysis/02_compare_visium_slideseq_cerebellum.md) -
  comparison of the Visium and Slide-seqV2 on different cerebellum
  samples
- [03_simulations.md](https://github.com/lulizou/spASE/blob/master/analysis/03_simulations.md) -
  simulation results
- [04_overall_results_summary_figures.md](https://github.com/lulizou/spASE/blob/master/analysis/04_overall_results_summary_figures.md) -
  main summary figures and tables from running spASE
- [05_celltype_xchr_figures.md](https://github.com/lulizou/spASE/blob/master/analysis/05_celltype_xchr_figures.md) -
  main figure comparing cell type-driven and within cell type spatial
  ASE
- [06_visium_slideseq_mixtures_overdispersion.md](https://github.com/lulizou/spASE/blob/master/analysis/06_visium_slideseq_mixtures_overdispersion.md) -
  supplemental figures showing count distributions of perfect binomial
  sampling, spatial transcriptomics, and Smart-seq3 data
- [07_129_CAST_analysis.md](https://github.com/lulizou/spASE/blob/master/analysis/07_129_CAST_analysis.md) -
  supplemental figures on comparing CAST and 129 transcripts
- [08_xchr_inactivation.md](https://github.com/lulizou/spASE/blob/master/analysis/08_xchr_inactivation.md) -
  supplemental figures showing X-chromosome genes for each sample

## Citation

The (old) pre-print is
[here](https://www.biorxiv.org/content/10.1101/2021.12.01.470861v1.full).
