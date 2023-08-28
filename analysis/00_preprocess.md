Preprocessing
================

- <a href="#sample-table" id="toc-sample-table">Sample table</a>
- <a href="#alignment" id="toc-alignment">Alignment</a>
- <a href="#rctd" id="toc-rctd">RCTD</a>
- <a href="#overall-ase-spase" id="toc-overall-ase-spase">Overall ASE:
  spASE</a>
- <a href="#cell-type-specific-analyses-c-side-and-spase"
  id="toc-cell-type-specific-analyses-c-side-and-spase">Cell type-specific
  analyses: C-SIDE and spASE</a>

This document details the pre-processing steps run before the rest of
the analysis scripts, including pointing to which scripts were run to
perform each step.

# Sample table

See code in .qmd document for LaTeX table.

# Alignment

For both Slide-seq and Visium data, we built custom bowtie2 indices for
the pooled transcriptome of the 129 and CAST mice.

Alignment was conducted with e.g. the command:

    bowtie2 -x bowtie2_index_129xCAST \
            -k 100 \
            -p 32 \
            --very-sensitive \
            -U ./tagged2.fastq |
            samtools view -bS - > ./tagged_bwt2_129_CAST.bam

BAM files were processed with a [custom Python
script](https://github.com/lulizou/spASE/blob/master/scripts/processBowtie2.py)
to get uniquely mapped reads based on number of mismatches, as well as
the number of reads uniquely mapped but unable to be assigned to one
allele (used as input for cell type assignment).

# RCTD

RCTD was run on each sample individually using the scripts (`>`
indicates the name of the output file produced from that script; note
that the scripts are not explicitly run with the pipe notation.):

- `run_rctd_hippo1.R > rctd_hippo_1.rds`
- `run_rctd_hippo2.R > rctd_hippo_2.rds`
- `run_rctd_hippo3.R > rctd_hippo_3.rds`
- `run_rctd_cere3.R > rctd_cere_3.rds`
- `run_rctd_cere4_visium.R > rctd_cere_4_visium.rds`
- `run_rctd_mix5_visium.R > rctd_mix_5_visium.rds`

See \`run_rctd.sbatch\`\` for SLURM job submission resources.

# Overall ASE: spASE

1.  Overall maternal/paternal bias - we assume
    $\text{logit}(p_{j}) = \beta_{0,j}$, i.e. the mean maternal
    probability does not change based on cell type or spatial location.

- `run_spase_hippo1_overall_bias.R > results_overall_bias_hippo_1.rds`
- `run_spase_hippo2_overall_bias.R > results_overall_bias_hippo_2.rds`
- `run_spase_hippo3_overall_bias.R > results_overall_bias_hippo_3.rds`
- `run_spase_cere3_overall_bias.R > results_overall_bias_cere_3.rds`
- `run_spase_cere4_visium_overall_bias.R > > results_overall_bias_cere_4_visium.rds`

See `run_spase_overall_bias.sbatch` for SLURM job submission resources.

2.  Overall spatial pattern (no cell type effect) - we assume
    $\text{logit}(p_{i,j}) = \beta_{0,j} + \sum_{\ell=1}^L x_{i,\ell}\beta_{\ell,j}$,
    where $x_{i,\ell}$ are degrees of freedom $L$ thin plate spline
    basis functions evaluated at spots $i$.

- `run_spase_hippo1_overall_spatial.R > results_overall_spatial_hippo_1.rds`
- `run_spase_hippo2_overall_spatial.R > results_overall_spatial_hippo_2.rds`
- `run_spase_hippo3_overall_spatial.R > results_overall_spatial_hippo_3.rds`
- `run_spase_cere3_overall_spatial.R > results_overall_spatial_cere_3.rds`
- `run_spase_cere4_visium_overall_spatial.R > results_overall_spatial_cere_4.rds`

See `run_spase_overall_spatial.sbatch` for SLURM job submission
resources.

# Cell type-specific analyses: C-SIDE and spASE

1.  Within cell type maternal/paternal bias - we assume
    $p_{i,j} = \sum_{k=1}^K \alpha_{i,j,k} \text{expit}(\beta_{0,k,j})$.
2.  Within cell type spatial pattern - we assume
    $p_{i,j} = \sum_{k=1}^K \alpha_{i,j,k} \, \text{expit} \left(\beta_{0,k,j} + \sum_{\ell=1}^L x_{i,\ell}\beta_{\ell,k,j}\right)$,
    where $x_{i,\ell}$ are degrees of freedom $L$ thin plate spline
    basis functions evaluated at spots $i$.

C-SIDE and spASE were run on each sample individually using the scripts:

- `run_spase_hippo1_celltype.R > results_celltype_hippo_1_df_5.rds`
- `run_spase_hippo2_celltype.R > results_celltype_hippo_2_df_5.rds`
- `run_spase_hippo3_celltype.R > results_celltype_hippo_3_df_5.rds`
- `run_spase_cere3_celltype.R > results_celltype_cere_3_df_5.rds`
- `run_spase_cere4_visium_celltype.R > results_celltype_cere_4_df_5.rds`

See `run_spase_celltype.sbatch` for example SLURM job submission
resources.
