#' Filtered and annotated contigs of TCR from mice
#'
#' The details of how these are generated are shown in the vignette mouse_tcell_qc
#' and are serialized to serve as an examples for other vignettes and documentation.
#' @format A data frame of 3399 contigs and 22 fields,
#'  all except 4 are originally defined in <https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/output/annotation#contig>
#'  The following fields were defined ex post facto.
#'
#'  1. `anno_file`: Path to original csv file
#'  2. `pop`: Mouse strain.
#'  3. `sample`: An artificial "replicate" from the original data defined by subsampling with replacement
#'  4. `celltype`: The putative cell type of the contig.
#'
"contigs_qc"

#' A preconstructed `ContigClusterDB` from the `contigs_qc` data
#'
#' Used in various examples.
#' @format
#' \code{ccdb_ex = ContigCellDB_10XVDJ(contigs_qc, contig_pk = c('pop',   'sample', 'barcode', 'contig_id'), cell_pk = c('pop',   'sample', 'barcode'))}
#'
"ccdb_ex"
