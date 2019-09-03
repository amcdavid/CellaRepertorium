#' Generate "pseudobulk" data from a `ContigCellDB`
#'
#' @param ccdb [ContigCellDB()]
#' @param class_keys `character` naming fields in `contig_tbl` that define unique classes of the repertoire
#' @param total_keys `character` naming fields that also define a class, and moreover will be conditioned upon when calculating the total
#' @param type one of "cell" or "umi"
#'
#' @return `tibble`
#' @export
#'
generate_pseudobulk = function(ccdb, class_keys, total_keys, type = c('cell', 'umi')){
    type = match.arg(type, c('cell', 'umi'))
    all_keys = unique(c(class_keys, total_keys))

    if(type == 'umi'){
        bulk <- ccdb$contig_tbl %>% group_by(!!!syms(all_keys)) %>% summarize(n_class = sum(umis))
    } else {
        bulk <- ccdb$contig_tbl %>% group_by(!!!syms(all_keys)) %>% summarize(n_class = dplyr::n_distinct(barcode))
    }

    bulk = bulk %>% ungroup() %>% tidyr::complete(.,!!!syms(total_keys), tidyr::nesting(!!!syms(class_keys)),fill=list(n_class=0))
    bulk <- bulk %>% group_by(!!!syms(total_keys)) %>% mutate(total = sum(n_class))
}

