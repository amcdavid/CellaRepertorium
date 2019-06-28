#' Filtering rows of tables in ContigCellDB object
#'
#' @param ccdb A ContigCellDB object
#' @param which_tbl name of the table needs to be filtered
#' @param ... conditions used to filter out rows
#'
#' @return ContigCellDB object with updated table
#' @export
#'
#' @examples
#' data(contigs_qc)
#' MIN_CDR3_AA <- 6
#'
#' cdb <- ContigCellDB_10XVDJ(contigs_qc, contig_pk = c('barcode', 'pop', 'sample', 'contig_id'), cell_pk = c('barcode', 'pop', 'sample'))
#' new_cdb <- filter_cdb(cdb,full_length, productive == 'True', high_confidence, chain != 'Multi', str_length(cdr3) > MIN_CDR3_AA)
filter_cdb <- function(ccdb, ..., which_tbl='contig_tbl'){
  tbl <- slot(ccdb,which_tbl)
  tbl <- dplyr::filter(.data=tbl,!!!rlang::quos(...))
  slot(ccdb,which_tbl) <- tbl
  return(ccdb)
}

#' Make changes (ceate new or update existed) on columns of tables in ContigCellDB object
#'
#' @param ccdb A ContigCellDB object
#' @param ... name and value pair of column that will be updated
#' @param which_tbl name of the table needs to update columns
#'
#' @return ContigCellDB object with updated table
#' @export
#'
#' @examples
#' data(contigs_qc)
#' MIN_CDR3_AA <- 6
#' cdb <- ContigCellDB_10XVDJ(contigs_qc, contig_pk = c('barcode', 'pop', 'sample', 'contig_id'), cell_pk = c('barcode', 'pop', 'sample'))
#' new_cdb <- mutate_cdb(cdb, new_col = 1)
mutate_cdb <- function(ccdb, ..., which_tbl='contig_tbl'){
  tbl <- slot(ccdb,which_tbl)
  tbl <- tbl %>% dplyr::mutate(!!!rlang::quos(...))
  slot(ccdb,which_tbl) <- tbl
  return(ccdb)
}


#' Implement cell QC
#'
#' @param contig_tbl contig tbl
#' @param contig_keys names of the slots in the tbl
#' @param plot logic variable. If TRUE, then plot the figure.
#'
#' @return a tbl or data frame
#' @export
#'
#' @examples
cellqc <- function(contig_tbl, contig_keys, plot = TRUE){
  total_umi <- contig_tbl %>% group_by(!!!syms(contig_keys)) %>% summarize(total_umi = sum(umis)) %>% rename(is_T = 'celltype == "T_ab')
  if (plot){
    ggplot(filter(total_umi, high_confidence, is_T), aes(color = factor(is_cell), x = total_umi, group = interaction(is_cell, sample, pop))) + stat_ecdf() + coord_cartesian(xlim = c(0, 10)) + ylab('Fraction of barcodes') + theme_minimal() + scale_color_discrete('10X called cell?')
  }
  return(total_umi)
}