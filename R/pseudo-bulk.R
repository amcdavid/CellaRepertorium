#' Generate "pseudobulk" data from a `ContigCellDB`
#'
#' Tabulate contigs with a unique combination of `class_keys` per `total_keys`.
#' For instance, `total_keys` might be a sample identifier, and `class_keys` might
#' be the V- and J- gene identities.  The idea is that this might mimic the data
#' generated in a bulk experiment.
#'
#' This function is currently rather 10x-specific, in that it is assumed that columns `barcode` and
#' `umis` exist.
#' @param ccdb [ContigCellDB()]
#' @param class_keys `character` naming fields in `contig_tbl` that define unique classes of the repertoire
#' @param total_keys `character` naming fields to be conditioned upon when calculating the total.
#' @param type one of "cell" or "umi"
#'
#' @return `tibble`
#' @export
#' @examples
#' data(ccdb_ex)
#' ccdb_ex = cluster_germline(ccdb_ex)
#' pseudo = generate_pseudobulk(ccdb_ex, c('v_gene', 'j_gene', 'chain'), c('pop', 'sample'))
generate_pseudobulk = function(ccdb, class_keys, total_keys, type = c('cell', 'umi')){
    type = match.arg(type, c('cell', 'umi'))
    all_keys = unique(c(class_keys, total_keys))
    total_keys = setdiff(total_keys, class_keys)
    if(!all(c('umis', 'barcode') %in% names(ccdb$contig_tbl))) stop("Expecting columns `barcode` and `umis` in `contig_tbl`.")
    contig_grp = ccdb$contig_tbl %>% group_by(!!!syms(all_keys))
    if(type == 'umi'){
        bulk = contig_grp %>% summarize(n_class = sum(.data$umis))
    } else {
        bulk = contig_grp %>% summarize(n_class = dplyr::n_distinct(.data$barcode))
    }

    bulk = ungroup(bulk)
    bulk = tidyr::complete(bulk,!!!syms(total_keys), tidyr::nesting(!!!syms(class_keys)),fill=list(n_class=0))
    bulk = group_by(bulk, !!!syms(total_keys)) %>% mutate(total = sum(.data$n_class))
}

