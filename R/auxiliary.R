#' Guess the cell type of a contig from the chain ID
#'
#' This function is likely dependent on annotations from 10X and
#' may change or break
#' as their pipeline changes.
#' @param chain `character` which will be parsed to try to infer celltype
#'
#' @return contig table with `celltype` column
#' @seealso [crosstab_by_celltype()]
#' @examples
#' data(ccdb_ex)
#' table(guess_celltype(ccdb_ex$contig_tbl$chain))
#' @export
guess_celltype = function(chain) {
      celltype = case_when(
        chain %in% c('TRA', 'TRB') ~ "T_ab",
        chain %in% c('TRD', 'TRG') ~ 'T_gd',
        chain == 'Multi' ~ 'Multi',
        chain %in% c('IGH', 'IGK', 'IGL') ~ 'B',
        TRUE ~ 'Others')
  celltype
}

#' Count contig UMIs by celltype
#'
#' @param ccdb A ContigCellDB object
#'
#' @return a table, keyed by `cell_pk` counting UMIs per celltype
#' @export
#' @seealso [guess_celltype()]
#' @examples
#' data(ccdb_ex)
#' nrow(ccdb_ex$cell_tbl)
#' total_umi = crosstab_by_celltype(ccdb_ex)
#' nrow(total_umi)
crosstab_by_celltype = function(ccdb) {
    check_contig_names(ccdb, c('chain', 'umis'))
    if('celltype__' %in% names(ccdb$cell_pk)) stop("Rename or remove `celltype__` from `cell_pk`.")
    ccdb$contig_tbl = ccdb$contig_tbl %>%
      dplyr::mutate(celltype__ = guess_celltype(.data$chain))

    # group by cell_keys
    cell_keys = union(ccdb$cell_pk, 'celltype__')
    total_umi = ccdb$contig_tbl %>% group_by(!!!syms(cell_keys)) %>%
      summarize(total_umi = sum(.data$umis))
    total_umi = left_join_warn(total_umi, ccdb$cell_tbl, by = ccdb$cell_pk, join = right_join) #in case cell_tbl has celltype__ column not as PK
    total_umi = tidyr::spread(total_umi, 'celltype__', 'total_umi', fill = 0)

    return(total_umi)
}
