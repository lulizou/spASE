# spASE
R package for detecting 2D spatial patterns of allele-specific expression in spatial
transcriptomics data. Provides functions for estimating, plotting, and
testing for significant ASE patterns.

## Installation

```r
devtools::install_github("lulizou/spASE")
```

## Details

The main function for estimating 2D smooth functions of ASE spatial data
is `spase`, which requires:
* matrix (rows = genes, columns = bead or spot barcodes) of UMI counts
  from allele 1
* matrix of UMI counts from allele 2 (same dimensions as first matrix)
* covariates data frame which contain the following columns (in order):
  1.  bead or spot barcodes
  2.  x1 coordinates of that bead or spot
  3.  x2 coordinates
  4.  (optional) factor covariate, such as cell type

Then the same information as well as the result of `spase` can be used in
`plotSpase` to generate smoothed allele probability functions, cross-sections
of those functions, and 2D z-score plots.

Code to reproduce figures and results can be found in the vignettes folder.
