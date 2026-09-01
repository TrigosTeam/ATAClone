# ATAClone

ATAClone uses single-cell ATAC-seq to identify tumour "clones" that share a copy number
profile, then estimates absolute copy number for each clone jointly.

## Installation

ATAClone is an R package and is not on CRAN. Install it from GitHub:

```r
# install.packages("remotes")
remotes::install_github("TrigosTeam/ATAClone")
```

Several dependencies come from Bioconductor and may need to be installed first:

```r
# install.packages("BiocManager")
BiocManager::install(c(
  "ComplexHeatmap", "GenomicRanges", "MatrixGenerics", "rtracklayer",
  "scran", "transformGamPoi", "glmGamPoi", "scDblFinder", "AUCell"
))
```

Other dependencies (`data.table`, `ggplot2`, `igraph`, `irlba`, `MASS`, `Matrix`,
`reshape2`, `uwot`, and `Seurat` for cell cycle gene sets) are available from CRAN.

## Getting started

The [multiome vignette](https://github.com/TrigosTeam/ATAClone/blob/main/inst/vignette/ataclone_multiome_vignette.Rmd)
demonstrates all key steps in an ATAClone analysis on a 10X multiome sample including: quality control and filtering, normalisation, clone identification,
doublet removal, and absolute copy number estimation.

## Input data

ATAClone works from matrices of genomic bin counts (bins × cells). The vignette uses two per
sample: one built from all ATAC-seq fragments, and one built from only fragments overlapping
stably accessible regions, which removes most of the technical variation in transposition
efficiency between droplets.

A normal reference dataset is needed for copy number estimation. This can be an external
non-tumour sample, or a non-tumour cluster from the same experiment when calling absolute
copy number. Matched RNA counts from a multiome assay are optional but help with quality
control, doublet detection, and cell cycle scoring.

Example datasets are installed with the package under `extdata/`.

## Contributing

Issues and pull requests are welcome at
<https://github.com/TrigosTeam/ATAClone/issues>.

## Citation

If you use ATAClone in your work, please cite:

> Cain L, Trigos AS. ATAClone: Cancer Clone Identification and Copy Number Estimation from
> Single-cell ATAC-seq. *bioRxiv* 2026.03.11.710984 (2026).
> doi: [10.64898/2026.03.11.710984](https://doi.org/10.64898/2026.03.11.710984)

```bibtex
@article{cain2026ataclone,
  title   = {ATAClone: Cancer Clone Identification and Copy Number Estimation from Single-cell ATAC-seq},
  author  = {Cain, Lachlan and Trigos, Anna S},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.03.11.710984},
  url     = {https://www.biorxiv.org/content/10.64898/2026.03.11.710984v2}
}
```
