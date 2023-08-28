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

The main input to spASE includes: - maternal counts matrix (rows are
genes, columns are spots or cells) - paternal counts matrix - spatial
coordinates.

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

To reproduce the results in the manuscript, see the `analysis` folder.
The data is available in the Broad single cell portal at
[SCP1692](https://singlecell.broadinstitute.org/single_cell/study/SCP1692/detection-of-allele-specific-expression-in-spatial-transcriptomics-with-spase).

## Citation

The (old) pre-print is
[here](https://www.biorxiv.org/content/10.1101/2021.12.01.470861v1.full).
