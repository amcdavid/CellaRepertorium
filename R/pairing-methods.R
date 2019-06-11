#' For each cell, return a single, canonical chain-cluster
#'
#' In single cell data, multiple chains (heavy-light or alpha-beta) are expected.  In some cases, there could be more than two (eg multiple alpha alleles for T cells).
#' This picks a cluster id for each cell based on the overall prevalence of cluster ids over all cells in `tbl`.
#' If order = 1 then the canonical chain-cluster will be the most prevalent, and if order = 2, it will be the 2nd most prevalent, and so on.  Ties are broken arbitrarily (possibly by lexicographic order of `cluster_idx`).
#' @param tbl `data.frame` containing columns specified in `cell_identifiers`, `cluster_idx` and optionally `chain_identifiers`
#' @param cell_identifiers `character` vector specifying columns in `tbl` that identify a cell
#' @param cluster_idx `character` specifying the column in `tbl` that identifies a clsuter
#' @param order return the 1st, 2nd, 3rd, etc, most common chain-cluster
#'
#' @return `data.frame` with columns from `cell_identifiers` and a single `cluster_idx` for each cell
#' @export
canonicalize_by_prevalence = function(tbl, cell_identifiers = 'barcode', cluster_idx = 'cluster_idx', order = 1){
    prevalence = tbl %>% group_by(!!sym(cluster_idx)) %>% summarize(prev = dplyr::n())
    tbl_prevalence = left_join(tbl, prevalence, by = cluster_idx)
    tbl_order = tbl_prevalence %>% group_by(!!!syms(cell_identifiers)) %>% summarize(!!cluster_idx := dplyr::nth(!!sym(cluster_idx), -order, prev))
    ungroup(tbl_order)
}




#' @param sort_factors `character` vector naming columns in `tbl` to sorted on, within  `cell_identifier`. Sorted by first element first, then ties broken by subsequent elements.  Sorted in decreasing order for each element.
#' @param chain_levels an optional `character` vector providing the sort order of the `chain` column in `tbl`.  Set to length zero to disable.
#' @export
#' @describeIn canonicalize_by_prevalence return a canonical contig by chain type, with TRB/IGH returned first. By default, ties are broken by umis and reads.
canonicalize_by_chain = function(tbl,  cell_identifiers = 'barcode', sort_factors = c('chain', 'umis', 'reads'), cluster_idx = 'cluster_idx', order = 1, chain_levels = c('IGL', 'IGK', 'TRA', 'TRB', 'IGH')){
    #ochain = tbl[[sort_factors[1]]]
    if(length(chain_levels) > 0){
        if(!('chain' %in% names(tbl))) stop('`tbl` must contain a `chain` column to set `chain_levels`')
        tbl = tbl %>% mutate(chain = factor(chain, levels = chain_levels, ordered = TRUE))
    }
    arranging = purrr::map(sort_factors, ~ rlang::quo(desc(!!sym(.x))))
    tbl %>% group_by(!!!syms(cell_identifiers)) %>% dplyr::arrange(!!!arranging) %>% mutate(rank = seq_along(!!sym(sort_factors[[1]]))) %>% dplyr::filter(rank == order) %>% dplyr::select(-rank)

}

# Define canonical chain types per cell

#' Given a family of similar sequences, return a "representative"
#'
#' @param seqs character vector
#' @param medoid_idx optional index into seqs
#' @param warn_if_distinct Should a warning be emitted if there are distinct elements in seqs?
#'
#' If medoid_idx is supplied, the medoid sequence is returned, otherwise the longest
#' sequence is returned
#'
#' @return character vector
#' @export
#'
#' @examples
#' get_canonical_representative(c('apple', 'manzana', 'pomme'))
get_canonical_representative = function(seqs, medoid_idx, warn_if_distinct = FALSE){
    if(!missing(medoid_idx)){
        if(!is.integer(medoid_idx) || medoid_idx < 1 || medoid_idx > length(seqs)) stop("Illegal `medoid_idx`")
        ii = medoid_idx
    } else{
        len = str_length(seqs)
        ii = which.max(len)
    }
    rep = seqs[ii]

    if(warn_if_distinct && length(nomatch <- which(seqs != rep)) > 0) warning("At indices ", nomatch, " sequences did not match the representative")

    return(rep)
}

#' Generate a list of tables representing clusters paired in cells
#'
#' A contingency table of every combination of `cluster_idx` up to `table_order` is generated.
#' Combinations that are found in at least `min_expansion` number of cells are reported.  All cells that have these combinations are returned, as well as cells that only have `orphan_level` of matching `cluster_idx`.
#'
#' For example, if `table_order=2` and `min_expansion=2` then heavy/light or alpha/beta pairs found two or more times will be returned (as well as alpha-alpha pairs, etc, if those are present).
#' If `orphan_level=1` then all cells that share just a single chain with an expanded clone will be returned.
#'
#' The `cluster_idx.1_fct` and `cluster_idx.2_fct` fields in `cell_tbl`, `idx1_tbl`, `idx2_tbl` are cast to factors and ordered such that pairings will tend to occur along the diagonal when they are cross-tabulated.
#' This facilitates plotting.
#'
#' @section Caveats and warnings:
#'  The cell_idx -> cluster_idx map is generally one-to-many, and is resolved by `canonicalize_fun`.  For `table_order>1`, few collisions are expected as most cells will contain no more than 2 clusters.  Any collisions are resolved by returning the most prevalent cluster, across samples.  When two clusters are tied for most prevalent within a cell, the `cluster_idx` returned is arbitrary. Therefore, when `table_order=1`, it is strongly recommended to subset the `cluster_tbl` to just a single chain.
#'
#' @param cluster_tbl a table with all combinations of clusters in all cells
#' @param cell_identifiers character vector naming fields that key a cell
#' @param cluster_idx character naming a single field IDing the clusters
#' @param min_expansion the minimal number of times a pairing needs to occur for it to be reported
#' @param cell_tbl optional, ancillary table with additional cell features.  Must also be keyed by `cell_identifiers`
#' @param feature_tbl optional, ancillary table with additional cluster features.  Must also be keyed by `cluster_idx`
#' @param table_order Integer larger than 1. What order of cluster_idx will be paired, eg, order = 2 means that the most common and second most common cluster_idx will be sought for each cell
#' @param orphan_level Integer larger than 0 and less than or equal to `table_order`.  Given that at least `min_expansion` cells are found that have `table_order` chains identical, how many `cluster_idx` pairs will we match on to select other cells.  Example: `ophan_level=1` means that cells that share just a single chain with the
#' @param cluster_whitelist a table of "cluster_idx" that should always be reported.  In contrast to the `cluster_tbl`, here the clusters must be named "cluster_idx.1", "cluster_idx.2" (if order-2 pairs are being selected).
#' @param cluster_blacklist a table of "cluster_idx" that will never be reported.  Must be named as per `cluster_whitelist`.
#' @param canonicalize_fun a function with signature `canonicalize_fun(cluster_tbl, cell_identifiers, cluster_idx, order = i)` that for each `cell_identifier` returns a single contig that depends on the `order`.  For instance \link{canonicalize_by_prevalence} or \link{canonicalize_by_chain}.
#'
#' @return list of tables.  The `cell_tbl` is keyed by the `cell_identifiers`, with fields "cluster_idx.1", "cluster_idx.2", etc, IDing the contigs present in each cell. "cluster_idx.1_fct" and "cluster_idx.2_fct" cast these fields to factors and are reordered to maximize the number of pairs along the diagonal. The `idx1_tbl` and `idx2_tbl` report information (passed in about the `cluster_idx` by `feature_tbl`.)  The `cluster_pair_tbl` reports all pairings found of contigs, and the number of times observed.
#' @export
#'
#' @seealso canonicalize_by_prevalence, canonicalize_by_chain
#' @importFrom tibble as_data_frame
#' @importFrom dplyr bind_rows left_join ungroup summarize anti_join
#' @importFrom stringr str_length str_c
#' @importFrom rlang sym syms :=
#' @examples
#' library(dplyr)
#' cluster_tbl = data_frame(clust_idx = gl(3, 2), cell_idx = rep(1:3, times = 2))
#' # no pairs found twice
#' pt1 = pairing_tables(cluster_tbl, 'cell_idx', 'clust_idx', canonicalize_by_prevalence)
#' # all pairs found, found once.
#' pt2 = pairing_tables(cluster_tbl, 'cell_idx', 'clust_idx',
#'     canonicalize_by_prevalence, min_expansion = 1)
#' pt2$cell_tbl
#' cluster_tbl2 = bind_rows(cluster_tbl, cluster_tbl %>% mutate(cell_idx = rep(4:6, times = 2)))
#' #all pairs found twice
#' pt3 = pairing_tables(cluster_tbl2, 'cell_idx', 'clust_idx', canonicalize_by_prevalence, min_expansion = 1)
#' pt3$cell_tbl
#' # `canonicalize_by_chain` expects fields `umis`, `reads`
#' # to break ties,  wrap the function to change this
#' cluster_tbl3 = cluster_tbl2 %>%
#'     mutate(umis = 1, reads = 1, chain = rep(c('TRA', 'TRB'), times = 6))
#' pt4 = pairing_tables(cluster_tbl3, 'cell_idx', 'clust_idx',
#'     canonicalize_by_chain, min_expansion = 1, table_order = 2)
pairing_tables = function(cluster_tbl, cell_identifiers = 'barcode', cluster_idx = 'cluster_idx', canonicalize_fun = canonicalize_by_chain, table_order = 2, min_expansion = 2,  orphan_level = 1, cluster_whitelist = NULL, cluster_blacklist = NULL, cell_tbl = NULL, feature_tbl = NULL ){

    if(orphan_level > table_order) stop('`ophan_level` must be less than or equal to `table_order`')
    if(table_order < 1) stop('Table order must be at least 1')

    # get `table_order` most common clusters for each cell
    # forcibly rename cluster_idx -> "cluster_idx"
    bar_chain_tbls = purrr::map(seq_len(table_order), function(i) canonicalize_fun(cluster_tbl, cell_identifiers = cell_identifiers, cluster_idx = cluster_idx, order = i) %>% dplyr::select(!!!c(syms(cell_identifiers), rlang::quo(cluster_idx))) %>% dplyr::rename( !!paste0('cluster_idx.', i) := !!cluster_idx))
    # for each cell, what clusters are present
    oligo_cluster_pairs = purrr::reduce(bar_chain_tbls, left_join, by = cell_identifiers)

    # Set up cluster_ids for indexing into the pairing tables
    cluster_ids = str_c('cluster_idx.', seq_len(table_order))
    cluster_ids_to_select = cluster_ids[seq_len(orphan_level)]

    # In how many cells do each cluster pairing appear?
    cluster_pair_tbl = oligo_cluster_pairs %>% group_by(!!!syms(cluster_ids)) %>% summarize(n_clone_pairs = dplyr::n())
    # which clusters are expanded
    expanded_cluster = cluster_pair_tbl %>% dplyr::filter(n_clone_pairs >= min_expansion) %>% dplyr::filter_at(.vars = cluster_ids, .vars_predicate = all_vars(!is.na(.)))
    expanded_cluster = ungroup(expanded_cluster) %>% dplyr::select(!!!syms(cluster_ids_to_select), max_pairs = n_clone_pairs)
    if(!is.null(cluster_whitelist)){
        expanded_cluster = bind_rows(expanded_cluster, cluster_whitelist)
    }
    if(!is.null(cluster_blacklist)) expanded_cluster = anti_join(expanded_cluster, cluster_blacklist)

    # Could have duplicated cluster_ids after binding to the whitelist or from considering orphans
    expanded_cluster = expanded_cluster[!duplicated(expanded_cluster %>% dplyr::select(-max_pairs)),]
    expanded_c1 = oligo_cluster_pairs %>% dplyr::inner_join(expanded_cluster, by = cluster_ids_to_select)
    if(anyDuplicated(expanded_c1 %>% dplyr::select(!!!syms(cell_identifiers)))) stop("Ruhoh, duplicated cell identifiers, this is a bug!")

    # cross-tab pairings in order to figure out how to order the cluster_idx to put most common pairing on diagonal
    if(nrow(expanded_c1)>0 && table_order > 1){
        expanded_counts = reshape2::dcast(expanded_c1, cluster_idx.1 ~ cluster_idx.2, fun.aggregate = length, value.var = 'max_pairs')
        mat = expanded_counts[,-1]
        #rownames(mat) = rowid[['cluster_idx']]
        ro = hclust(dist(mat))$order
        co = hclust(dist(t(mat)))$order
    } else if (nrow(expanded_c1) == 0) {
        warning('No pairs found')
        return()
    } else{
        expanded_counts = reshape2::dcast(expanded_c1, cluster_idx.1 ~ ., fun.aggregate = length, value.var = 'max_pairs')
        ro = order(expanded_counts[[2]])
        co = NA_integer_
    }

    ci_class = class(cluster_tbl[[cluster_idx]])
    as_method = if(ci_class == 'factor') as.factor else function(x) as(x, ci_class)
    rowid = data_frame(cluster_idx = expanded_counts$cluster_idx.1 %>% as_method, plot_order = ro)
    colid = suppressWarnings(data_frame(cluster_idx = colnames(expanded_counts)[-1] %>% as_method, plot_order = co))
    rowid[['cluster_idx.1_fct']] = factor(rowid[['cluster_idx']], levels = rowid[['cluster_idx']][ro])
    colid[['cluster_idx.2_fct']] = factor(colid[['cluster_idx']], levels = colid[['cluster_idx']][co])


    # also fix levels of this tbl
    expanded_c1 = expanded_c1 %>% mutate(cluster_idx.1_fct = factor(cluster_idx.1, levels = levels(rowid[['cluster_idx.1_fct']])))
    if(table_order>1) expanded_c1 = expanded_c1 %>% mutate(cluster_idx.2_fct = factor(cluster_idx.2, levels = levels(colid[['cluster_idx.2_fct']])))

    if(!is.null(cell_tbl)){
        if(anyDuplicated(cell_tbl %>% select(!!!syms(cell_identifiers)))) stop('`cell_tbl` must not have duplicate `cell_identifiers`.')
        expanded_c1 = left_join(expanded_c1, cell_tbl, by = cell_identifiers)
    }

    idx1_tbl = rowid %>% dplyr::rename(!!cluster_idx := cluster_idx)
    idx2_tbl = colid %>% dplyr::rename(!!cluster_idx := cluster_idx)
    if(!is.null(feature_tbl)){
        if(anyDuplicated(feature_tbl %>% select(!!sym(cluster_idx)))) stop('`feature_tbl` must not have duplicate `cluster_idx`.')
        idx1_tbl = left_join(idx1_tbl, feature_tbl, by = cluster_idx)
        idx2_tbl = left_join(idx2_tbl, feature_tbl, by = cluster_idx)
    }

    list(cell_tbl = expanded_c1, idx1_tbl = idx1_tbl, idx2_tbl = idx2_tbl, cluster_pair_tbl = cluster_pair_tbl)

}

plot_pairing = function(pairing_list, color_labels_by){
    pl = pairing_list
    pairs_plt = ggplot(pairing_list$cell_tbl, aes(x = cluster_idx.1_fct, y = cluster_idx.2_fct, color = sample, shape = pop)) + geom_jitter(width = .3, height = .3)

    ylab = data_frame(!!color_labels_by :=  ggplot_build(pairs_plt)$layout$panel_params[[1]]$y.label) %>% left_join(feature_tbl) %>% mutate(class_color = ifelse(is.na(class_color), '#E41A1C', class_color))

    xlab = data_frame(!!color_labels_by :=  ggplot_build(pairs_plt)$layout$panel_params[[1]]$x.label) %>% left_join(feature_tbl) %>% mutate(class_color = ifelse(is.na(class_color), '#E41A1C', class_color))

    pairs_plt = pairs_plt + theme(axis.text.x = element_text(angle = 90, color = xlab$class_color, size = 8), axis.text.y = element_text(color = ylab$class_color, size = 8))


}


ig_chain_recode = function(tbl){
    pairing = case_when(tbl$IGH>0 & (tbl$IGK>0 | tbl$IGL>0) ~ 'paired',
                        tbl$IGH>0 ~ 'heavy',
                        (tbl$IGK>0 | tbl$IGL>0) ~ 'light')
    canonical = case_when(tbl$IGH<2 & (tbl$IGK==2 | tbl$IGL==2) ~ 'double-light',
                          tbl$IGH<2 & ((tbl$IGK + tbl$IGL)>1) ~ 'multi-light',
                          tbl$IGH<2 & (tbl$IGK + tbl$IGL)<2 ~ 'classical',
                          tbl$IGH>1 ~ 'multi-heavy',
                          TRUE ~ 'other')
    dplyr::bind_cols(tbl, tibble(pairing, canonical))
}

tcr_chain_recode = function(tbl){
    pairing = case_when(tbl$TRA>0 & tbl$TRB>0 ~ 'paired',
                        tbl$TRB>0 ~ 'beta',
                        tbl$TRA>0 ~ 'alpha')
    canonical = case_when(tbl$TRB==2 ~ 'double-beta',
                          tbl$TRA==2  ~ 'double-alpha',
                          (tbl$TRB + tbl$TRA) > 1 ~ 'other',
                          TRUE ~ 'classical')
    dplyr::bind_cols(tbl, tibble(pairing, canonical))
}

#' Categorize the pairing present in a cell
#'
#' For each cell (defined by `ccdb$cell_pk`) count the number of each level of `chain_key` occurs, and cross tabulate.
#' Also for each cell, paste together all values `chain_key`.
#' Return a tibble, keyed by cells that includes the counts of the chains, the `raw_chain_type` and any additional output from running `chain_recode_fun`.
#' @param ccdb `ContigCellDB`
#' @param chain_key `character` naming the field in the `contig_tbl` identifying chain
#' @param chain_recode_fun a function that operates on the output of this function that further reduces the chain combinations to some other summary.  Set to 'guess' to apply functions that may work for 10X data or `NULL` to skip.  See `CellaRepertorium:::tcr_chain_recode` for an example.
#'
#' @return a `tibble` keyed by cells.
#' @export
#'
#' @examples
#' data(ccdb_ex)
#' enumerate_pairing(ccdb_ex)
#' enumerate_pairing(ccdb_ex, chain_recode_fun = 'guess')
enumerate_pairing = function(ccdb, chain_key = 'chain', chain_recode_fun = NULL){
    if(!is.null(chain_recode_fun) && chain_recode_fun == 'guess'){
        top_chain = names(sort(table(ccdb$contig_tbl[[chain_key]]), decreasing = TRUE))[1]
        if(top_chain %in% c('TRA', 'TRB')){
            chain_recode_fun = tcr_chain_recode
        } else if(top_chain %in% c('IGH', 'IGK', 'IGL')){
            chain_recode_fun = ig_chain_recode
        }
    } else if(is.null(chain_recode_fun)){
        chain_recode_fun = function(x) x
    }
    if(!is.function(chain_recode_fun)) stop("`chain_recode_fun` must be a function, NULL, or 'guess'")

    chain_keys = union(chain_key, ccdb$cell_pk)
    chain_count = ccdb$contig_tbl %>% group_by(!!!syms(chain_keys)) %>% summarize(n_chains = n()) %>% spread(chain_key, 'n_chains', fill = 0)
    chain_type = ccdb$contig_tbl %>% group_by(!!!syms(ccdb$cell_pk)) %>% summarize(raw_chain_type = paste(sort(!!sym(chain_key)), collapse = '_'))
    chain_summary = left_join(chain_type, chain_count, by = ccdb$cell_pk) %>% ungroup()
    chain_recode_fun(chain_summary)
}
