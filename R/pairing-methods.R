#' Rank contigs, per cell, by experiment-wide prevalence of `cluster_pk`, which is added as the `prevalence` field
#'
#' @return `ContigCellDB` with modified `contig_tbl`
#' @inheritParams canonicalize_cell
#' @export
#' @examples
#' data(ccdb_ex)
#' ccdb_ex = cluster_germline(ccdb_ex)
#' rank_prev = rank_prevalence_ccdb(ccdb_ex)
#' rank_prev$contig_tbl
#' rank_chain = rank_chain_ccdb(ccdb_ex)
#' rank_chain$contig_tbl
rank_prevalence_ccdb = function(ccdb, contig_filter_args = TRUE,
                              tie_break_keys = c('umis', 'reads')){
    ccdb$contig_tbl = ccdb$contig_tbl %>% group_by(!!!syms(ccdb$cluster_pk)) %>% mutate(prevalence = dplyr::n())
    ranked_contigs = .rank_contigs_per_cell(ccdb, tie_break_keys = c('prevalence', tie_break_keys), rlang::enexpr(contig_filter_args))
   ccdb$contig_tbl =  ranked_contigs
   ccdb
}




#' @param chain_levels an optional `character` vector providing the sort order of the `chain` column in `tbl`.  If set to length zero, then the the ordering will be alphabetical
#' @param chain_key `character` naming the field in `contig_tbl` to be sorted on.
#' @describeIn rank_prevalence_ccdb return a canonical contig by chain type, with TRB/IGH returned first. By default, ties are broken by umis and reads.
#' @export
rank_chain_ccdb = function(ccdb,  contig_filter_args = TRUE,
                           tie_break_keys = c('umis', 'reads'),
                           chain_key = 'chain',
                           contig_fields = tie_break_keys, chain_levels = c('IGL', 'IGK', 'TRA', 'TRB', 'IGH')){
    check_contig_names(ccdb, chain_key)
    if(length(chain_levels) > 0){
        ccdb$contig_tbl[[chain_key]] = factor(ccdb$contig_tbl[[chain_key]],  levels = chain_levels, ordered = TRUE)
    }
    ranked_contigs = .rank_contigs_per_cell(ccdb, tie_break_keys = c(chain_key, tie_break_keys), rlang::enexpr(contig_filter_args))
    ccdb$contig_tbl = ranked_contigs
    ccdb
}


check_contig_names = function(ccdb, req_names) {
  if (length(missing_contig <- setdiff(req_names, names(ccdb$contig_tbl))) > 0) stop('`contig_tbl` is missing fields, ', paste(missing_contig, collapse = ', '), '.')
    return(TRUE)
}

.rank_contigs_per_cell <- function(ccdb, tie_break_keys, contig_filter_args) {
    tbl = ccdb$contig_tbl
    check_contig_names(ccdb, tie_break_keys)
    # Filter with expressions in contig_filter_args
    #filter_arg = rlang::enexpr(contig_filter_args)
    ft = filter(.data = tbl, !!contig_filter_args)
    # setup quosures to arrange the data
    arranging = purrr::map(tie_break_keys, ~ rlang::quo(desc(!!sym(.x))))
    # take first row of each cell
    ft2 = ft %>% group_by(!!!syms(ccdb$cell_pk)) %>% dplyr::arrange(!!!arranging)
    ranked_contigs = ft2 %>% dplyr::mutate(n_grp = dplyr::n(), grp_rank = seq_along(.data$n_grp))
    ungroup(ranked_contigs)
}

#' Find a canonical contig to represent a cell
#'
#' Using filtering in `contig_filter_args` and sorting in `tie_break_keys` and `order` find a
#' single, canonical contig to represent each cell
#' Fields in `contig_fields` will be copied over to the `cell_tbl`.
#' @inheritParams canonicalize_cluster
#'
#' @return [ContigCellDB()] with some number of clusters/contigs/cells but with "canonical" values copied into `cell_tbl`
#' @export
#' @seealso [canonicalize_cluster()]
#' @examples
#' # Report beta chain with highest umi-count, breaking ties with reads
#' data(ccdb_ex)
#' beta = canonicalize_cell(ccdb_ex, chain == 'TRB',
#' tie_break_keys = c('umis', 'reads'),
#' contig_fields = c('umis', 'reads', 'chain', 'v_gene', 'd_gene', 'j_gene'))
#' head(beta$cell_tbl)
#'
#' # Stable: only adds fields to `cell_tbl`
#' stopifnot(dplyr::all_equal(beta$cell_tbl[ccdb_ex$cell_pk],
#' ccdb_ex$cell_tbl[ccdb_ex$cell_pk], ignore_row_order = TRUE))
#'
#' #Report cdr3 with highest UMI count, but only when > 5 UMIs support it
#' umi5 = canonicalize_cell(ccdb_ex, umis > 5,
#' tie_break_keys = c('umis', 'reads'), contig_fields = c('umis', 'cdr3'))
#' stopifnot(all(umi5$cell_tbl$umis > 5, na.rm = TRUE))
canonicalize_cell = function(ccdb, contig_filter_args = TRUE,
                             tie_break_keys = c('umis', 'reads'),
                             contig_fields = tie_break_keys, order = 1, overwrite = TRUE){
    check_contig_names(ccdb, contig_fields)
    ranked_contigs = .rank_contigs_per_cell(ccdb, tie_break_keys, rlang::enexpr(contig_filter_args))
    ft2 = ranked_contigs[ranked_contigs$grp_rank==order,,drop = FALSE]
    cell_tbl = ccdb$cell_tbl
    # join with cell tbl (so same number of cells)
    ccdb$cell_tbl = right_join_warn(ft2[unique(c(contig_fields, ccdb$cell_pk))], cell_tbl, by = ccdb$cell_pk, overwrite = overwrite)
    ccdb
}

# This should return a ccdb or a subclass?
#' Generate a list of tables representing clusters paired in cells
#'
#' A contingency table of every combination of `cluster_idx` up to `table_order`
#' is generated. Combinations that are found in at least `min_expansion` number
#' of cells are reported.  All cells that have these combinations are returned,
#' as well as cells that only have `orphan_level` of matching `cluster_idx`.
#'
#' For example, if `table_order=2` and `min_expansion=2` then heavy/light or
#' alpha/beta pairs found two or more times will be returned
#' (as well as alpha-alpha pairs, etc, if those are present).
#' If `orphan_level=1` then all cells that share just a single chain with an
#' expanded clone will be returned.
#'
#' The `cluster_idx.1_fct` and `cluster_idx.2_fct` fields in `cell_tbl`,
#' `idx1_tbl`, `idx2_tbl` are cast to factors and ordered such that pairings will
#' tend to occur along the diagonal when they are cross-tabulated.
#' This facilitates plotting.
#'
#' @param ccdb `ContigCellDB`
#' @param ranking_key field in `ccdb$contig_tbl` giving the ranking of each contig per cell.  Probably generated by a call to [rank_prevalence_ccdb()] or [rank_chain_ccdb()].
#' @param min_expansion the minimal number of times a pairing needs to occur for
#' it to be reported
#' @param cluster_keys optional `character` naming additional columns in
#' `ccdb$cluster_tbl` to be reported in the pairing
#' @param table_order Integer larger than 1. What order of cluster_idx will be
#' paired, eg, order = 2 means that the first and second highest ranked contigs will be sought and paired in each cell
#' @param orphan_level Integer in interval \[1, `table_order`\].  Given that at least `min_expansion` cells are found that have `table_order` chains identical, how many `cluster_idx` pairs will we match on to select other cells.  Example: `ophan_level=1` means that cells that share just a single chain with an expanded pair will be reported.
#' @param cluster_whitelist a table of pairings or clusters that should always be reported.  Here the clusters must be named "cluster_idx.1", "cluster_idx.2" (if order-2 pairs are being selected) rather than with `ccdb$cluster_pk``
#' @param cluster_blacklist a table of pairings or clusters that will never be reported.  Must be named as per `cluster_whitelist`.
#'
#' @return list of tables.  The `cell_tbl` is keyed by the `cell_identifiers`, with fields "cluster_idx.1", "cluster_idx.2", etc, IDing the contigs present in each cell. "cluster_idx.1_fct" and "cluster_idx.2_fct" cast these fields to factors and are reordered to maximize the number of pairs along the diagonal. The `idx1_tbl` and `idx2_tbl` report information (passed in about the `cluster_idx` by `feature_tbl`.)  The `cluster_pair_tbl` reports all pairings found of contigs, and the number of times observed.
#' @export
#'
#' @seealso [rank_prevalence_ccdb()]
#' @importFrom tibble as_data_frame tibble
#' @importFrom dplyr bind_rows left_join ungroup summarize anti_join
#' @importFrom stringr str_length str_c
#' @importFrom rlang sym syms :=
#' @examples
#' library(dplyr)
#' tbl = tibble(clust_idx = gl(3, 2), cell_idx = rep(1:3, times = 2), contig_idx = 1:6)
#' ccdb = ContigCellDB(tbl, contig_pk = c('cell_idx', 'contig_idx'),
#' cell_pk = 'cell_idx', cluster_pk = 'clust_idx')
#' # add `grp_rank` to ccdb$contig_tbl indicating how frequent a cluster is
#' ccdb = rank_prevalence_ccdb(ccdb, tie_break_keys = character())
#' # using `grp_rank` to determine pairing
#' # no pairs found twice
#' pt1 = pairing_tables(ccdb)
#' # all pairs found, found once.
#' pt2 = pairing_tables(ccdb, min_expansion = 1)
#' pt2$cell_tbl
#' tbl2 = bind_rows(tbl, tbl %>% mutate(cell_idx = rep(4:6, times = 2)))
#' ccdb2 = ContigCellDB(tbl2, contig_pk = c('cell_idx', 'contig_idx'), cell_pk = 'cell_idx',
#' cluster_pk = 'clust_idx') %>% rank_prevalence_ccdb(tie_break_keys = character())
#' #all pairs found twice
#' pt3 = pairing_tables(ccdb2, min_expansion = 1)
#' pt3$cell_tbl
#' ccdb2$contig_tbl = ccdb2$contig_tbl %>%
#'     mutate(umis = 1, reads = 1, chain = rep(c('TRA', 'TRB'), times = 6))
#' ccdb2 = rank_chain_ccdb(ccdb2, tie_break_keys = character())
#' pt4 = pairing_tables(ccdb2, min_expansion = 1, table_order = 2)
pairing_tables = function(ccdb, ranking_key = 'grp_rank', table_order = 2, min_expansion = 2,  orphan_level = 1, cluster_keys = character(), cluster_whitelist = NULL, cluster_blacklist = NULL){

    if(orphan_level > table_order) stop('`ophan_level` must be less than or equal to `table_order`')
    if(table_order < 1) stop('Table order must be at least 1')

    # get `table_order` most common clusters for each cell
    # forcibly rename cluster_idx -> "cluster_idx"
    contig_tbl = ccdb$contig_tbl
    contig_tbl = contig_tbl[contig_tbl[[ranking_key]]<=table_order,]
    table_order = max(contig_tbl[[ranking_key]])
    cell_identifiers = ccdb$cell_pk
    cluster_idx = ccdb$cluster_pk
    cell_tbl = ccdb$cell_tbl
    cell_contig_tab = tidyr::pivot_wider(contig_tbl[c(cell_identifiers, cluster_idx, ranking_key)], names_from = !!ranking_key, values_from = cluster_idx)

    # Set up cluster_ids for indexing into the pairing tables
    cluster_ids =  str_c('cluster_idx.', seq_len(table_order))
    names(cell_contig_tab)[(ncol(cell_contig_tab)-table_order+1):ncol(cell_contig_tab)] = cluster_ids
    cluster_ids_to_select = cluster_ids[seq_len(orphan_level)]

    # In how many cells do each cluster pairing appear?
    cluster_pair_tbl = cell_contig_tab %>% group_by(!!!syms(cluster_ids)) %>% summarize(n_clone_pairs = dplyr::n())
    # which clusters are expanded,
    expanded_cluster = dplyr::filter(cluster_pair_tbl, .data$n_clone_pairs >= min_expansion)
    # Must have both cluster_ids non NA (otherwise not a pairing of the required order)
    expanded_cluster = dplyr::filter_at(expanded_cluster, .vars = cluster_ids, .vars_predicate = dplyr::all_vars(!is.na(.)))
    expanded_cluster = ungroup(expanded_cluster) %>% dplyr::select(!!!syms(cluster_ids_to_select), max_pairs = .data$n_clone_pairs)
    if(!is.null(cluster_whitelist)){
        expanded_cluster = bind_rows(expanded_cluster, cluster_whitelist)
    }
    if(!is.null(cluster_blacklist)) expanded_cluster = anti_join(expanded_cluster, cluster_blacklist)

    # Could have duplicated cluster_ids after binding to the whitelist or from considering orphans
    expanded_cluster = expanded_cluster[!duplicated(expanded_cluster %>% dplyr::select(-.data$max_pairs)),]
    expanded_c1 = cell_contig_tab %>% dplyr::inner_join(expanded_cluster, by = cluster_ids_to_select)
    if(anyDuplicated(expanded_c1[cell_identifiers])) stop("Ruhoh, duplicated cell identifiers, this is a bug!")

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

    ci_class = class(contig_tbl[[cluster_idx]])
    as_method = if(ci_class == 'factor') as.factor else function(x) as(x, ci_class)
    rowid = tibble(cluster_idx = expanded_counts$cluster_idx.1 %>% as_method, plot_order = ro)
    colid = suppressWarnings(tibble(cluster_idx = colnames(expanded_counts)[-1] %>% as_method, plot_order = co))
    rowid[['cluster_idx.1_fct']] = factor(rowid[['cluster_idx']], levels = rowid[['cluster_idx']][ro])
    colid[['cluster_idx.2_fct']] = factor(colid[['cluster_idx']], levels = colid[['cluster_idx']][co])


    # also fix levels of this tbl
    expanded_c1 = expanded_c1 %>% mutate(cluster_idx.1_fct = factor(.data$cluster_idx.1, levels = levels(rowid[['cluster_idx.1_fct']])))
    if(table_order>1) expanded_c1 = expanded_c1 %>% mutate(cluster_idx.2_fct = factor(.data$cluster_idx.2, levels = levels(colid[['cluster_idx.2_fct']])))

    if(!is.null(cell_tbl)){
        if(anyDuplicated(cell_tbl %>% select(!!!syms(cell_identifiers)))) stop('`cell_tbl` must not have duplicate `cell_identifiers`.')
        expanded_c1 = left_join(expanded_c1, cell_tbl, by = cell_identifiers)
    }

    idx1_tbl = rowid %>% dplyr::rename(!!cluster_idx := cluster_idx)
    idx2_tbl = colid %>% dplyr::rename(!!cluster_idx := cluster_idx)
    feature_tbl = ccdb$cluster_tbl[unique(c(cluster_keys, cluster_idx))]
    idx1_tbl = left_join_warn(idx1_tbl, feature_tbl, by = cluster_idx)
    idx2_tbl = left_join_warn(idx2_tbl, feature_tbl, by = cluster_idx)

    list(cell_tbl = expanded_c1, idx1_tbl = idx1_tbl, idx2_tbl = idx2_tbl, cluster_pair_tbl = cluster_pair_tbl)

}


#' @export
#' @describeIn enumerate_pairing Recode a table with IG chains
#' @importFrom dplyr case_when
#' @importFrom tibble tibble
ig_chain_recode = function(tbl){
  pairing = case_when(tbl$IGH>0 & (tbl$IGK>0 | tbl$IGL>0) ~ 'paired',
                      tbl$IGH>0 ~ 'heavy',
                      (tbl$IGK>0 | tbl$IGL>0) ~ 'light',
                      tbl$IGH==0 & tbl$IGK==0 & tbl$IGL ==0 ~ 'none')
  canonical = case_when(tbl$IGH<2 & (tbl$IGK + tbl$IGL)<2 ~ 'classical',
      tbl$IGH>2 & (tbl$IGK + tbl$IGL)<2 ~ 'multi-heavy',
      tbl$IGH<2 & (tbl$IGK + tbl$IGL) == 2 ~ 'double-light',
      tbl$IGH<2 & ((tbl$IGK + tbl$IGL)>2) ~ 'multi-light',
      tbl$IGH == 2 & (tbl$IGK + tbl$IGL)<2 ~ 'double-heavy',
      tbl$IGH == 0 & (tbl$IGK + tbl$IGL) == 0 ~ 'none',
      TRUE ~ 'other')
  dplyr::bind_cols(tbl, tibble(pairing, canonical))
}

#' @export
#' @describeIn enumerate_pairing Recode a table with TCR chains
#' @param tbl output from enumerate_pairing containing TRA/TRB or IGH/IHK/IHL columns
tcr_chain_recode = function(tbl){
    pairing = case_when(tbl$TRA>0 & tbl$TRB>0 ~ 'paired',
                        tbl$TRB>0 ~ 'beta',
                        tbl$TRA>0 ~ 'alpha',
                        tbl$TRA == 0 & tbl$TRB == 0 ~ 'none')
    canonical = case_when(tbl$TRB==2 ~ 'double-beta',
                          tbl$TRA==2  ~ 'double-alpha',
                          (tbl$TRB + tbl$TRA) > 2 ~ 'other',
                          tbl$TRA == 0 & tbl$TRB == 0 ~ 'none',
                          TRUE ~ 'classical')
    dplyr::bind_cols(tbl, tibble(pairing, canonical))
}



#' Generate a 2d cross tab using arbitrary numbers of columns as factors
#'
#' As many rows as unique combs of x_fields
#' As many columns as unique combs of y_fields
#' No NA.
#' @param tbl `data.frame`
#' @param x_fields `character` fields in `tbl`
#' @param y_fields `character` fields in `tbl`
#'
#' @return `tibble`
#' @export
#'
#' @examples
#' cross_tab_tbl(mtcars, c('cyl', 'gear'), 'carb')
cross_tab_tbl = function(tbl, x_fields, y_fields){
    x_key = unique(tbl[x_fields]) %>% mutate(x_key__ = factor(seq_len(nrow(.))))
    y_key = unique(tbl[y_fields]) %>% mutate(y_key__ = factor(seq_len(nrow(.))))
    tbl = left_join(tbl, x_key) %>% left_join(y_key)
    cross_tb = unclass(table(tbl[['x_key__']], tbl[['y_key__']], exclude = NULL))
    cross_tbl = as_tibble(cross_tb)
    names(cross_tbl) = do.call(paste, c(y_key[y_fields], list(sep = '_')))
    cross_tbl[['x_key__']] = factor(seq_len(nrow(cross_tbl)))
    left_join(cross_tbl, x_key) %>% select(-x_key__)
}

#' Categorize the pairing present in a cell
#'
#' For each cell (defined by `ccdb$cell_pk`) count the number of each level of `chain_key` occurs, and cross tabulate.
#' Also for each cell, paste together all values `chain_key`.
#' Return a tibble, keyed by cells that includes the counts of the chains, the `raw_chain_type` and any additional output from running `chain_recode_fun`.
#' @param ccdb `ContigCellDB`
#' @param chain_key `character` naming the field in the `contig_tbl` identifying chain
#' @param chain_recode_fun a function that operates on the output of this function that further reduces the chain combinations to some other summary.  Set to 'guess' to apply functions that may work for 10X data or `NULL` to skip.  See `CellaRepertorium::tcr_chain_recode` for an example.
#'
#' @return a `tibble` keyed by cells.
#' @export
#'
#' @examples
#' data(ccdb_ex)
#' enumerate_pairing(ccdb_ex)
#' enumerate_pairing(ccdb_ex, chain_recode_fun = 'guess')
enumerate_pairing = function(ccdb, chain_key = 'chain', chain_recode_fun = NULL){
    if(!is.null(chain_recode_fun) && !is.function(chain_recode_fun) && chain_recode_fun == 'guess'){
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
    tbl = left_join_warn(ccdb$cell_tbl[ccdb$cell_pk], ccdb$contig_tbl, by = ccdb$cell_pk)
    tbl[[chain_key]] = forcats::fct_explicit_na(tbl[[chain_key]], na_level = 'none')
    chain_count = tbl %>% group_by(!!!syms(chain_keys)) %>% summarize(n_chains = dplyr::n())
    chain_table = tidyr::spread(chain_count, chain_key, 'n_chains', fill = 0)
    chain_table = chain_table[setdiff(names(chain_table), 'none')]
    chain_type = tbl %>% group_by(!!!syms(ccdb$cell_pk)) %>% summarize(raw_chain_type = paste(sort(!!sym(chain_key)), collapse = '_'))
    chain_summary = left_join(chain_table, chain_type, by = ccdb$cell_pk) %>% ungroup()
    chain_summary = left_join_warn(ccdb$cell_tbl[ccdb$cell_pk], chain_summary, by = ccdb$cell_pk)
    chain_recode_fun(chain_summary)
}
