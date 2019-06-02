#' Filtered and annotated contigs of TCR from mice
#'
#' The details of how these are generated are shown in the vignette mouse_tcell_qc
#' and are serialied to serve as an examples for other vignettes and documentation.
#' @format A data frame of 3399 contigs and 22 fields,
#'  all except 4 are originally defined in https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/output/annotation#contig
#'  The following fields were defined ex post facto
#'  \describe{
#'    \item{anno_file}{Path to original csv file}
#'    \item{pop}{Mouse strain}
#'    \item{sample}{An artificial "replicate" from the original data defined by subsampling with replacement}
#'    \item{celltype}{The putative cell type of the contig}
#'  }
"contigs_qc"
