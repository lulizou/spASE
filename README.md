# spASE
Detecting 2D spatial patterns of allele-specific expression in spatial
transcriptomics data. Provides functions for estimating, plotting, and
testing for significant ASE patterns.

## Installation

```r
devtools::install_github("lulizou/spASE")
```

## Details

The main function for estimating 2D smooth functions of ASE spatial data
is `spase`, which requires:
\begin{itemize}
  \item matrix (rows = genes, columns = bead or spot barcodes) of UMI counts
  from allele 1
  \item matrix of UMI counts from allele 2 (same dimensions as first matrix)
  \item covariates which contain the following columns (in order):
  \begin{enumerate}
    \item bead or spot barcodes
    \item x1 coordinates of that bead or spot
    \item x2 coordinates
  \end{enumerate}
\end{itemize}

Code to reproduce figures and results can be found in the vignettes folder.
