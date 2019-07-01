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


#' make the cross table by celltype
#'
#' @param ccdb A ContigCellDB object
#'
#' @return a cross table
#' @export
#'
#' @examples
#' cdb = ContigCellDB(all_anno,contig_pk = c('barcode','pop','sample','contig_id'),cell_pk = c('barcode','pop','sample'))
#' total_umi <- crosstab_by_celltype(cdb)
crosstab_by_celltype <- function(ccdb){
  # add celltype column
  ccdb$contig_tbl <- ccdb$contig_tbl %>% dplyr::mutate(celltype = case_when(chain %in% c('TRA', 'TRB') ~ "T_ab", chain %in% c('TRD', 'TRG') ~ 'T_gd', chain == 'Multi' ~ 'Multi', chain %in% c('IGH','IGK', 'IGL') ~ 'B', TRUE ~ NA_character_))
  
  # group by cell_keys
  cell_keys <- union(ccdb$cell_pk,'celltype')
  total_umi <- ccdb$contig_tbl %>% group_by(!!!syms(cell_keys)) %>% summarize(total_umi = sum(umis)) %>% tidyr::spread(celltype, 'total_umi', fill = 0)
  
  
  return(total_umi)
}